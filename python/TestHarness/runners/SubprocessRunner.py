#* This file is part of the MOOSE framework
#* https://www.mooseframework.org
#*
#* All rights reserved, see COPYRIGHT for full restrictions
#* https://github.com/idaholab/moose/blob/master/COPYRIGHT
#*
#* Licensed under LGPL 2.1, please see LICENSE for details
#* https://www.gnu.org/licenses/lgpl-2.1.html

import os, platform, subprocess, shlex, time
from tempfile import SpooledTemporaryFile
from signal import SIGTERM
from TestHarness.runners.Runner import Runner
from TestHarness import util

class SubprocessRunner(Runner):
    """
    Runner that spawns a local subprocess.
    """
    def __init__(self, job, options):
        Runner.__init__(self, job, options)

        # The output file handler
        self.outfile = None
        # The error file handler
        self.errfile = None
        # The underlying subprocess
        self.process = None
        # The joined output
        self.output = None

    def spawn(self, timer):
        tester = self.job.getTester()
        use_shell = tester.specs["use_shell"]
        cmd = tester.getCommand(self.options)
        tester.setCommandRan(cmd)

        # Split command into list of args to be passed to Popen
        if not use_shell:
            cmd = shlex.split(cmd)

        self.process = None
        self.outfile = SpooledTemporaryFile(max_size=1000000) # 1M character buffer
        self.errfile = SpooledTemporaryFile(max_size=100000)  # 100K character buffer

        process_args = [cmd]
        process_kwargs = {'stdout': self.outfile,
                            'stderr': self.errfile,
                            'close_fds': False,
                            'shell': use_shell,
                            'cwd': tester.getTestDir()}
        # On Windows, there is an issue with path translation when the command is passed in
        # as a list.
        if platform.system() == "Windows":
            process_kwargs['creationflags'] = subprocess.CREATE_NEW_PROCESS_GROUP
        else:
            process_kwargs['preexec_fn'] = os.setsid

        # Set this for OpenMPI so that we don't clobber state
        if tester.hasOpenMPI():
            popen_env = os.environ.copy()
            popen_env['OMPI_MCA_orte_tmpdir_base'] = tester.getTempDirectory().name
            process_kwargs['env'] = popen_env

        try:
            self.process = subprocess.Popen(*process_args, **process_kwargs)
        except Exception as e:
            raise Exception('Error in launching a new task') from e

        timer.start()

    def wait(self, timer):
        self.process.wait()

        timer.stop()

        self.exit_code = self.process.poll()
        self.outfile.flush()
        self.errfile.flush()

        # store the contents of output, and close the file
        self.output = util.readOutput(self.outfile, self.errfile, self.job.getTester())
        self.outfile.close()
        self.errfile.close()

    def kill(self):
        if self.process is not None:
            try:
                if platform.system() == "Windows":
                    from distutils import spawn
                    if spawn.find_executable("taskkill"):
                        subprocess.call(['taskkill', '/F', '/T', '/PID', str(self.process.pid)])
                    else:
                        self.process.terminate()
                else:
                    pgid = os.getpgid(self.process.pid)
                    os.killpg(pgid, SIGTERM)
            except OSError: # Process already terminated
                pass

    def getOutput(self):
        return self.output

    def sendSignal(self, signal):
        # process.poll() returns the process's exit code if it has completed,
        # and None if it is still running. This acts as a safety precaution
        # against an infinite loop; this will always close.
        while self.process.poll() is None:

            # tell() gives the current position in the file. If it is greater
            # than zero, the binary has started running and writing output. If
            # the output is blank, the moose_test binary hasn't actually started
            # doing anything yet. If so, sleep briefly and check again.
            if not self.outfile.tell():
                time.sleep(0.05)

            # If the output isn't blank, then we finally send the signal and exit the loop
            else:
                os.kill(self.process.pid, signal)
                break
