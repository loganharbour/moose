#* This file is part of the MOOSE framework
#* https://www.mooseframework.org
#*
#* All rights reserved, see COPYRIGHT for full restrictions
#* https://github.com/idaholab/moose/blob/master/COPYRIGHT
#*
#* Licensed under LGPL 2.1, please see LICENSE for details
#* https://www.gnu.org/licenses/lgpl-2.1.html

import itertools, re, os, time, threading
from timeit import default_timer as clock
from TestHarness.StatusSystem import StatusSystem
from TestHarness.FileChecker import FileChecker
from TestHarness.runners.Runner import Runner
from TestHarness import util
from tempfile import TemporaryDirectory
import traceback

class Timer(object):
    """
    A helper class for testers to track the time it takes to run.

    Every call to the start method must be followed by a call to stop.
    """
    def __init__(self):
        self.starts = []
        self.ends = []
    def start(self):
        """ starts the timer clock """
        self.starts.append(clock())
    def stop(self):
        """ stop/pauses the timer clock """
        self.ends.append(clock())
    def cumulativeDur(self):
        """ returns the total/cumulative time taken by the timer """
        diffs = [end - start for start, end in zip(self.starts, self.ends)]
        return sum(diffs)
    def averageDur(self):
        return self.cumulativeDur() / len(self.starts)
    def nRuns(self):
        return len(self.starts)
    def reset(self):
        self.starts = []
        self.ends = []

class Job(object):
    """
    The Job class is a simple container for the tester and its associated output file object, the DAG,
    the process object, the exit codes, and the start and end times.
    """
    # Iterator for producing a unique Job ID
    id_iter = itertools.count()

    def __init__(self, tester, job_dag, options):
        self.id = next(self.id_iter)
        self.options = options
        self.__j_lock = threading.Lock()
        self.__tester = tester
        self.specs = tester.specs
        self.__job_dag = job_dag
        self.timer = Timer()
        self.__start_time = clock()
        self.__end_time = None
        self.__previous_time = None
        self.__joined_out = ''
        self.report_timer = None
        self.__slots = None
        self.__meta_data = {}

        # Create a fileChecker object to be able to call filecheck methods
        self.fileChecker = FileChecker(self.options.input_file_name)

        # List of files modified by this job.
        self.modifiedFiles = []

        # Set of all jobs that this job races with.
        self.racePartners = set()

        # Alternate text we want to print as part of our status instead of the
        # pre-formatted status text (LAUNCHED (pbs) instead of FINISHED for example)
        self.__job_message = ''

        ### Enumerate the job statuses we want to use
        self.job_status = StatusSystem()

        self.hold = self.job_status.hold
        self.queued = self.job_status.queued
        self.running = self.job_status.running
        self.skip = self.job_status.skip
        self.error = self.job_status.error
        self.timeout = self.job_status.timeout
        self.finished = self.job_status.finished

        self.__finished_statuses = [self.skip,
                                    self.error,
                                    self.timeout,
                                    self.finished]

        self.__pending_statuses = [self.hold,
                                   self.queued,
                                   self.running]

        # Initialize jobs with a holding status
        self.setStatus(self.hold)

        # Whether or not we should forcefully report the status of this Job
        # the next time report statuses
        self.force_report_status = False

        # The object that'll actually do the run
        self._runner = None

        # Any additional output produced by the Job (not from the Tester or Runner)
        self.output = ''

        self.cached_output = None

        # A temp directory for this Job, if requested
        self.tmp_dir = None

    def __del__(self):
        # Do any cleaning that we can (removes the temp dir for now if it exists)
        self.cleanup()

    def getID(self):
        """Returns the unique ID for the job"""
        return self.id

    def setRunner(self, runner: Runner):
        """Sets the underlying Runner object that will run the command"""
        self._runner = runner

    def getLock(self):
        """ Get the lock associated with this job """
        return self.__j_lock

    def getTempDirectory(self):
        """
        Gets a shared temp directory that will be cleaned up for this Tester
        """
        if self.tmp_dir is None:
            self.tmp_dir = TemporaryDirectory(prefix='tester_')
        return self.tmp_dir

    def cleanup(self):
        """
        Entry point for doing any cleaning if necessary.

        Currently just cleans up the temp directory
        """
        if self.tmp_dir is not None:
            # Don't let this fail
            try:
                self.tmp_dir.cleanup()
            except:
                pass
            self.tmp_dir = None

    def getUpstreams(self):
        """ Return a list of all the jobs that needed to be completed before this job """
        dag = self.getDAG()
        original_dag = dag.getOriginalDAG()
        return dag.predecessors(self, original_dag)

    def getDownstreams(self):
        """ Return a list of all the jobs that need this job completed """
        dag = self.getDAG()
        original_dag = dag.getOriginalDAG()
        return dag.all_downstreams(self, original_dag)

    def getDAG(self):
        """ Return the DAG associated with this tester """
        return self.__job_dag.getDAG()

    def getTester(self):
        """ Return the tester object """
        return self.__tester

    def getSpecs(self):
        """ Return tester spec params """
        return self.getTester().specs

    def getTestName(self):
        """ Wrapper method to return the testers test name """
        return self.__tester.getTestName()

    def getTestNameShort(self):
        """ Return the shorthand Test name """
        return self.__tester.getTestNameShort()

    def getPrereqs(self):
        """ Wrapper method to return the testers prereqs """
        return self.__tester.getPrereqs()

    def getTestDir(self):
        """ Wrapper method to return the testers working directory """
        return self.__tester.getTestDir()

    def addCaveats(self, kwargs):
        """ Wrapper method for setting caveats """
        return self.__tester.addCaveats(kwargs)

    def getCaveats(self):
        """ Wrapper method for getting caveats """
        return self.__tester.getCaveats()

    def clearCaveats(self):
        """ Wrapper method for clearing caveats """
        return self.__tester.clearCaveats()

    def getCommand(self):
        """ Wrapper method for returing command """
        return self.__tester.getCommand(self.options)

    def getCommandRan(self):
        """ Wrapper method for returing command ran """
        return self.__tester.getCommandRan()

    def getRunnable(self):
        """ Wrapper method to return getRunnable """
        return self.__tester.getRunnable(self.options)

    def getOutputFiles(self, options):
        """ Wrapper method to return getOutputFiles """
        return self.__tester.getOutputFiles(options)

    def getMaxTime(self):
        """ Wrapper method to return getMaxTime """
        return self.__tester.getMaxTime()

    def getInputFile(self):
        """ Wrapper method to return input filename """
        return self.__tester.getInputFile()

    def getInputFileContents(self):
        """ Wrapper method to return input file contents """
        return self.__tester.getInputFileContents()

    def getUniqueIdentifier(self):
        """ A unique identifier for this job object """
        return os.path.join(self.getTestDir(), self.getTestName())

    def getUniquePrereqs(self):
        """ Return a list of prereqs with what should be their unique identifier """
        prereqs = self.getPrereqs()
        unique_prereqs = []
        for prereq in prereqs:
            unique_prereqs.append(os.path.join(self.getTestDir(), prereq))
        return unique_prereqs

    def addMetaData(self, **kwargs):
        """ Allow derived methods to store additional data which ends up in the data storage file """
        for key, value in kwargs.items():
            self.__meta_data[key] = value

    def getMetaData(self):
        """ return data stored by addMetaData """
        return self.__meta_data

    def getSlots(self):
        """ Return the number of slots this job consumes """
        if self.__slots == None:
            return self.setSlots(self.__tester.getSlots(self.options))
        return self.__slots

    def setSlots(self, slots):
        """ Set the number of processors this job consumes """
        self.__slots = int(slots)
        return self.__slots

    def canParallel(self):
        """ Call DAG and determine if this group can run in parallel """
        return self.__job_dag.canParallel()

    def run(self):
        """
        A blocking method to handle the exit status of the process object while keeping track of the
        time the process was active. When the process exits, read the output and close the file.
        """
        tester = self.__tester

        # Helper for trying and catching
        def try_catch(do, exception_name):
            try:
                do()
            except:
                self.cleanup()
                self.setStatus(self.error, f'{exception_name} EXCEPTION')
                self.output += util.outputHeader('Python exception encountered')
                self.output += traceback.format_exc()
                return False
            return True

        # Do not execute app, but still run the tester
        # This is truly awful and I really hate that it got put in here,
        # please remove it if you can.
        if not tester.shouldExecute():
            run_tester = lambda: tester.run(self.options, 0, '')
            try_catch(run_tester, 'TESTER RUN')
            return

        if self.options.pedantic_checks and self.canParallel():
            # Before the job does anything, get the times files below it were last modified
            self.fileChecker.get_all_files(self, self.fileChecker.getOriginalTimes())
            self.addCaveats('pedantic check')
            time.sleep(1)

        tester.prepare(self.options)

        # Verify that the working directory is available right before we execute
        if not os.path.exists(tester.getTestDir()):
            self.setStatus(self.error, 'WORKING DIRECTORY NOT FOUND')
            return
        # Getting the command can also cause a failure, so try that
        tester.getCommand(self.options)
        if tester.isError():
            return

        self.timer.reset()

        self.__start_time = clock()

        # Spawn the process
        spawn = lambda: self._runner.spawn(self.timer)
        if not try_catch(spawn, 'RUNNER SPAWN'):
            return

        # Entry point for testers to do other things
        post_spawn = lambda: tester.postSpawn(self._runner)
        if not try_catch(post_spawn, 'TESTER POST SPAWN'):
            return

        # And wait for it to complete
        wait = lambda: self._runner.wait(self.timer)
        if not try_catch(wait, 'RUNNER WAIT'):
            return

        self.__start_time = self.timer.starts[0]
        self.__end_time = self.timer.ends[-1]

        # Job error occurred, which means the Runner didn't complete
        # so don't process anything else
        if self.isError():
            self.cleanup()
            return

        # And do finalize (really just cleans up output)
        self._runner.finalize()

        if self.options.pedantic_checks and self.canParallel():
            # Check if the files we checked on earlier were modified.
            self.fileChecker.get_all_files(self, self.fileChecker.getNewTimes())
            self.modifiedFiles = self.fileChecker.check_changes(self.fileChecker.getOriginalTimes(), self.fileChecker.getNewTimes())

        # Allow derived proccessResults to process the output and set a failing status (if it failed)
        runner_output = self._runner.getOutput()
        exit_code = self._runner.getExitCode()
        run_tester = lambda: tester.run(self.options, exit_code, runner_output)
        try_catch(run_tester, 'TESTER RUN')

        # Run cleanup now that we're done
        self.cleanup()

    def killProcess(self):
        """ Kill remaining process that may be running """
        if self._runner:
            try:
                self._runner.kill()
            except:
                pass
        self.cleanup()

    def getStartTime(self):
        """ Return the time the process started """
        return self.__start_time

    def getEndTime(self):
        """ Return the time the process exited """
        return self.__end_time

    def getOutput(self):
        """ Return the combined contents of output """
        # Cached output is used when reading from a results file,
        # when we don't run anything and just populate results
        if self.cached_output:
            return self.cached_output

        # Concatenate output in order of Runner, Tester, Job
        output = ''
        object_outputs = [self.getRunner().getOutput() if self.getRunner() else '',
                          self.getTester().getOutput() if self.getTester() else '',
                          self.output]
        for object_output in object_outputs:
            if object_output:
                # Append an extra line if we're missing one
                if output and output[-1] != '\n':
                    output += '\n'
                output += object_output

        return output

    def getRunner(self):
        """ Gets the Runner that actually runs the command """
        return self._runner

    def getOutputFile(self):
        """ Return the output file path """
        if ((self.options.ok_files
             or self.options.fail_files
             or self.options.sep_files)
             and (self.isPass() or self.isFail())):
            (status, message, color, exit_code, sort_value) = self.getJointStatus()
            output_dir = self.options.output_dir if self.options.output_dir else self.getTestDir()
            output_file = os.path.join(output_dir,
                                       '.'.join([os.path.basename(self.getTestDir()),
                                                 self.getTestNameShort().replace(os.sep, '.'),
                                                 status,
                                                 'txt']))
            return os.path.join(output_dir, output_file)

    def appendOutput(self, output):
        self.output += output

    def getActiveTime(self):
        """ Return active time """
        m = re.search(r"Active time=(\S+)", self.getOutput())
        if m != None:
            return float(m.group(1))

    def getSolveTime(self):
        """ Return solve time """
        m = re.search(r"solve().*", self.getOutput())
        if m != None:
            return m.group().split()[5]

    def setPreviousTime(self, t):
        """
        Allow an arbitrary time to be set. This is used by the QueueManager
        to set the time as recorded by a previous TestHarness instance.
        """
        self.__previous_time = t

    def getTiming(self):
        """ Return active time if available, if not return a comparison of start and end time """
        if self.getActiveTime():
            return self.getActiveTime()
        elif self.getEndTime() and self.getStartTime():
            return self.timer.cumulativeDur()
        elif self.getStartTime() and self.isRunning():
            # If the test is still running, return current run time instead
            return max(0.0, clock() - self.getStartTime())
        elif self.__previous_time:
            return self.__previous_time
        else:
            return 0.0

    def getStatus(self):
        return self.job_status.getStatus()

    def setStatus(self, status, message=''):
        if self.isFinished():
            return self.getStatus()
        # Try and preserve job messages
        if not self.__job_message and message:
            self.__job_message = message
        return self.job_status.setStatus(status)

    def createStatus(self):
        return self.job_status.createStatus()

    def previousTesterStatus(self, options, previous_storage=None):
        return self.__tester.previousTesterStatus(options, previous_storage)

    def getStatusMessage(self):
        return self.__job_message

    ### Boolean status comparisons based on current Job _and_ Tester status. All finalized status calls
    ### should now call job.isSomeStatus for the definitive answer.
    # the following are more job related...
    def isError(self):
        return self.getStatus() in self.job_status.getFailingStatuses()

    def isSkip(self):
        _status = self.getStatus()
        return (_status == self.finished and self.__tester.isSkip()) \
            or (_status == self.skip and self.isNoStatus()) \
            or (_status == self.skip and self.__tester.isSkip())

    def isHold(self):
        return self.getStatus() == self.hold
    def isQueued(self):
        _status = self.getStatus()
        return (_status == self.queued and self.isNoStatus()) \
            or (_status in self.__finished_statuses and self.__tester.isQueued())
    def isRunning(self):
        return self.getStatus() in self.job_status.getPendingStatuses()
    def isTimeout(self):
        return self.getStatus() == self.timeout
    def isPending(self):
        return self.getStatus() in self.__pending_statuses
    def isFinished(self):
        return self.getStatus() in self.__finished_statuses

    # the following more tester related...
    def isSilent(self):
        return self.__tester.isSilent()
    def isNoStatus(self):
        return self.__tester.isNoStatus()
    def isSilent(self):
        return self.__tester.isSilent() or (not self.options.report_skipped and self.isSkip())
    def isPass(self):
        return self.__tester.isPass()
    def isFail(self):
        return self.__tester.isFail() or self.isError()
    def isDiff(self):
        return self.__tester.isDiff()
    def isDeleted(self):
        return self.__tester.isDeleted()

    def getJointStatus(self):
        """
        Return the most accurate message possible, based on the statuses of both the Tester and Job.
        Most of the time, we want to return a Tester status. As this class will have the most detail.
        """
        # Job has failed, or tester has no status
        if self.isError() or self.isNoStatus():
            return (self.getStatus().status,
                    self.getStatusMessage(),
                    self.getStatus().color,
                    self.getStatus().code,
                    self.getStatus().sort_value)

        # Tester has a finished status of some sort
        else:
            return (self.__tester.getStatus().status,
                    self.__tester.getStatusMessage(),
                    self.__tester.getStatus().color,
                    self.__tester.getStatus().code,
                    self.__tester.getStatus().sort_value)
