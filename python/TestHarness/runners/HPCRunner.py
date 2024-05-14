#* This file is part of the MOOSE framework
#* https://www.mooseframework.org
#*
#* All rights reserved, see COPYRIGHT for full restrictions
#* https://github.com/idaholab/moose/blob/master/COPYRIGHT
#*
#* Licensed under LGPL 2.1, please see LICENSE for details
#* https://www.gnu.org/licenses/lgpl-2.1.html

from TestHarness.runners.Runner import Runner
import re, time, os, subprocess

class HPCRunner(Runner):
    """
    Base Runner to be used with HPC schedulers (PBS, slurm)
    """
    def __init__(self, job, options, run_hpc):
        super().__init__(job, options)

        # The RunHPC object
        self.run_hpc = run_hpc

        # Number of seconds to try to wait for the output
        # We don't want to wait forever for output because
        # if the job ended in an unexpected state, it might
        # not even be using the output and we don't want to
        # just hang forever
        self.wait_output_time = 120

    def spawn(self, timer):
        self.run_hpc.submitJob(self.job)
        timer.start()

    def wait(self, timer):
        # Need to import here to avoid cyclic includes
        from TestHarness.schedulers.RunHPC import RunHPC

        # Poll loop waiting for the job to be finished
        # This gets a structure that represents the job, and the
        # polling itself is only done on occasion within RunHPC
        while True:
            time.sleep(1)
            hpc_job = self.run_hpc.getHPCJob(self.job)

            # We're done
            if hpc_job.done:
                self.exit_code = hpc_job.exit_code
                break

        timer.stop()

        # The PBS output (stdout+stderr)
        output_file = self.run_hpc.getHPCJobOutputPath(self.job)

        # If the Job is already finished, something happened in PBS
        # so we have an invalid state for processing in the Tester
        if self.job.isFinished():
            self.exit_code = -1

            # If we have output, we should try to add it
            self.trySetOutput()

            # Don't bother looking for the rest of the output
            return

        tester = self.job.getTester()

        # We've actually ran something now and not just qsub, so update the
        # command to what was ran there
        tester.setCommandRan(hpc_job.command)

        # Determine the output files that we need to wait for to be complete
        wait_files = set([output_file])
        # Output files needed by the Tester, only if it says we should
        if tester.mustOutputExist():
            for file in tester.getOutputFiles(self.options):
                wait_files.add(os.path.join(tester.getTestDir(), file))
        # The files that we can read, but are incomplete (no terminator)
        incomplete_files = set()

        # Wait for all of the files to be available
        file_poll_interval = 0.25
        waited_time = 0
        while wait_files or incomplete_files:
            # Look for each file
            for file in wait_files.copy():
                if os.path.exists(file) and os.path.isfile(file):
                    wait_files.discard(file)
                    incomplete_files.add(file)

            # Check for file completeness
            for file in incomplete_files.copy():
                if self.fileIsReady(file):
                    # Store the output
                    if file == output_file:
                        self.trySetOutput(throw=True)
                    # Done with this file
                    incomplete_files.discard(file)

            # We've waited for files for too long
            if (wait_files or incomplete_files) and waited_time >= self.wait_output_time:
                self.job.setStatus(self.job.error, 'FILE TIMEOUT')
                if self.output is None:
                    self.trySetOutput()
                def print_files(files, type):
                    if files:
                        self.output += '#' * 80 + f'\n{type} output file(s)\n' + '#' * 80 + '\n'
                        for file in files:
                            self.output += file + '\n'
                        self.output += '\n'
                print_files(wait_files, 'Unavailable')
                print_files(incomplete_files, 'Incomplete')
                break

            waited_time += file_poll_interval
            time.sleep(file_poll_interval)

    def kill(self):
        self.run_hpc.killJob(self.job)

    def trySetOutput(self, throw=False):
        if self.output is None:
            self.output = ''

        output_file = self.run_hpc.getHPCJobOutputPath(self.job)
        if os.path.exists(output_file) and os.path.isfile(output_file):
            try:
                # If we're trying to read output, we can't truncate it
                if self.job.getTester().needFullOutput(self.options):
                    self.output = open(output_file, 'r').read()
                else:
                    self.output = self.readTruncated(output_file)

                # If we can parse the exit code here, do it. Sometimes PBS
                # will do screwy stuff with not capturing the actual exit code...
                find_exit_code = re.search('Completed TestHarness RunHPC test execution; exit code = (\d+)', self.output)
                if find_exit_code:
                    self.exit_code = int(find_exit_code.group(1))
            except:
                if throw:
                    raise
                pass

    def fileIsReady(self, file):
        """
        Checks if a file is ready for reading.

        In summary:
        - Check if the file exists
        - If the file exists, make sure that it has the terminator
          string (to know that we have the full file)
        - Remove the terminator string
        """
        # The file terminator check (to we have the up-to-date copy of the file)
        # is dependent on whether or not the file is a binary
        is_binary = self.isFileBinary(file)
        # If this returns None, it means that the "file" command couldn't determine
        # the file type, which may be the case if we have an incomplete file so
        # just continue and check on the next iteration
        if is_binary is None:
            return False

        ending_comment = self.run_hpc.getOutputEndingComment()

        # Binary file
        if is_binary:
            with open(file, "rb+") as file:
                # We'll be looking for this many characters
                len_comment = len(ending_comment)

                # Move to the end and figure out the position
                # back where our terminator should be
                file.seek(0, os.SEEK_END)
                pos = file.tell() - len_comment

                # File is shorter than our comment
                if pos < 0:
                    return False

                # Move to the position where our terminator _should_ be
                file.seek(pos)

                # We try here in the event that we're loading
                # an earlier part of the file and we can't decode
                try:
                    contents = file.read(len_comment).decode('utf-8')
                except:
                    return False

                # Terminator isn't there
                if contents != ending_comment:
                    return False

                # Remove the terminator
                file.seek(pos)
                file.truncate()

                return True
        # Text file
        else:
            line, pos = self.getLastLine(file)
            if ending_comment == line:
                with open(file, "r+", encoding="utf-8") as f:
                    f.seek(pos)
                    f.truncate()
                return True

        return False

    @staticmethod
    def isFileBinary(file):
        """
        Returns whether or not the given file is a binary file.

        If None, a failure was encountered when checking the file type.
        """
        try:
            call_file = subprocess.check_output(['file', '--mime-encoding', file], text=True)
        except:
            return None

        # Will return something like "<filename>: <encoding>",
        # where <encoding>=binary when the file is binary
        find_binary = re.search('binary$', call_file)
        return find_binary is not None

    @staticmethod
    def getLastLine(file):
        """
        Gets the last line of a text file and the position
        in the file at which that last line is.
        """
        with open(file, 'rb') as f:
            try:
                f.seek(-2, os.SEEK_END)
                while f.read(1) != b'\n':
                    f.seek(-2, os.SEEK_CUR)
            except OSError: # one line filecd
                f.seek(0)
            pos = f.tell()
            line = f.readline().decode('utf-8')
            return line, pos

    @staticmethod
    def readTruncated(file, start_lines=1000, end_lines=1000):
        """
        Reads a file and truncates it past a certain amount of lines.
        """
        with open(file, 'rb') as f:
            # Find the end position of the file so that we don't read past
            f.seek(0, os.SEEK_END)
            total_bytes = f.tell()

            # Read the set of lines
            f.seek(0)
            head_lines_read = 0
            head = ''
            while head_lines_read < start_lines and f.tell() < total_bytes:
                head += f.read(1).decode('utf-8')
                if len(head) > 1 and head[-1:] == '\n':
                    head_lines_read += 1

            # Keep the end of the head position so that we don't read
            # backwards past it for the tail
            head_pos = f.tell()

            # Seek to the end and start reading ending lines
            f.seek(0, os.SEEK_END)

            # Keep reading the ending lines until we've reached the max
            # number of lines we want or have reached the head output
            tail_lines_read = 0
            tail = []
            while tail_lines_read < end_lines and f.tell() > head_pos:
                # Read each character in the line until we reach
                # the beginning or a new line
                line = []
                while f.tell() > 1:
                    f.seek(-2, os.SEEK_CUR)
                    char = f.read(1).decode('utf-8')
                    if char == '\n' or f.tell() == 0:
                        break
                    line.append(char)

                # Append the new read line
                line.reverse()
                tail.append(''.join(line))
                tail_lines_read += 1

            # Whether or not we have truncated output
            # (have hit the location of the head output)
            truncated = f.tell() != head_pos

        # Form the combined output
        output = head
        if truncated:
            output += f'{"#" * 80}\nOUTPUT TRIMMED\n{"#" * 80}\n'
        if tail:
            tail.reverse()
            output += '\n'.join(tail)

        return output
