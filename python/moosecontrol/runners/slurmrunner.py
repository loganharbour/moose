# This file is part of the MOOSE framework
# https://mooseframework.inl.gov
#
# All rights reserved, see COPYRIGHT for full restrictions
# https://github.com/idaholab/moose/blob/master/COPYRIGHT
#
# Licensed under LGPL 2.1, please see LICENSE for details
# https://www.gnu.org/licenses/lgpl-2.1.html

"""Defines the SlurmRunner."""

import os
from logging import getLogger
from re import match
from time import sleep
from typing import Optional, Tuple
import subprocess

import paramiko

from .socketrunner import SocketRunner

logger = getLogger("SlurmRunner")


class SlurmRunner(SocketRunner):
    """
    Base interface for subprocess runners.

    Is an abstract base class and is used in
    SubprocessSocketRunner and SubprocessPortRunner.
    """

    import paramiko

    def __init__(self, slurm_host: str, slurm_script: str, *args, **kwargs):

        assert isinstance(slurm_host, str)
        assert isinstance(slurm_script, str)

        # The slurm host to connect to
        self._slurm_host: str = slurm_host

        # Path to the slurm submission script on the remote
        self._slurm_script: str = slurm_script

        self._ssh_client: Optional[paramiko.SSHClient] = None

        self._job_id: Optional[int] = None

        super().__init__(socket_path=SocketRunner.random_socket_path(), *args, **kwargs)

    @property
    def slurm_host(self) -> str:
        """Get the slurm host."""
        return self._slurm_host

    @property
    def slurm_script(self) -> str:
        """Get the path to the slurm script on the remote."""
        return self._slurm_script

    def get_ssh_client(self) -> paramiko.SSHClient:
        if self._ssh_client is None:
            key_filename = None
            try:
                ssh_config = os.path.expanduser("~/.ssh/config")
                config = paramiko.SSHConfig.from_path(ssh_config).lookup(
                    self.slurm_host
                )
                identityfile = config.get("identityfile")
                if identityfile is not None and len(identityfile) > 0:
                    key_filename = identityfile[-1]
            except:
                pass

            client = paramiko.SSHClient()
            client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            try:
                client.connect(self.slurm_host, key_filename=key_filename)
            except Exception as e:
                raise ConnectionError(
                    f"Failed to connect to host {self.slurm_host}"
                ) from e

            self._ssh_client = client

        assert isinstance(self._ssh_client, paramiko.SSHClient)
        return self._ssh_client

    def exec_ssh_command(self, command: str) -> Tuple[int, str]:
        """
        Execute the given command on the slurm host.

        Parameters
        ----------
        command : str
            The command to run.

        Returns
        -------
        int :
            The exit code from the command.
        str :
            The on-screen output from the command.

        """
        logger.debug(f"Sending SSH command {command}")
        client = self.get_ssh_client()
        try:
            _, stdout, stderr = client.exec_command(command)
        except Exception as e:
            raise RuntimeError(
                f'Failed to run SSH command "{command}" on {self.slurm_host}'
            ) from e
        exit_code = stdout.channel.recv_exit_status()
        result = "".join(stdout.readlines())
        if exit_code != 0:
            result += "".join(stderr.readlines())
        return exit_code, result.rstrip()

    @staticmethod
    def error_exit(message: str):
        raise SystemExit(f"ERROR: {message}")

    def get_job_state(self) -> Tuple[str, Optional[str]]:
        """
        Get the state of the running Slurm job.

        Should only be ran after the job is submitted.

        Returns
        -------
        str : The job state.

        """
        job_id = self._job_id
        assert isinstance(job_id, int)

        cmd = f"sacct -j {job_id} --parsable2 --noheader -o state,NodeList"
        exit_code, result = self.exec_ssh_command(cmd)

        if exit_code != 0:
            raise RuntimeError(
                f"Failed to get job status on {self.slurm_host}; result:"
                f"\n\n{result}"
            )

        # The result (after it is running) will have multiple lines,
        # of which the first one is the status for the main job. So, first,
        # split on lines and just capture the first line
        result = result.splitlines()[0]
        # After this, we'll either have:
        # "STATUS|None assigned" if the node isn't assigned yet, or
        # "STATUS|node1234,node1235" if the node is assigned
        result_split = result.splitlines()[0].split("|")
        assert len(result_split) == 2
        # First | delimited entry, state
        state = result_split[0]
        # Second delimited entry, node list (if any)
        node = (
            None
            if result_split[1].startswith("None")
            else result_split[1].split(",")[0]
        )

        return state, node

    def initialize(self, data: dict):
        slurm_host = self.slurm_host
        slurm_script = self.slurm_script
        socket_path = self.socket_path

        # Setup the SSH client
        self.get_ssh_client()

        # Make sure the slurm path exists
        if self.exec_ssh_command(f"test -f {slurm_script}")[0] != 0:
            raise FileNotFoundError(
                f"Slurm script {slurm_script} on {slurm_host} does not exist"
            )

        # Submit the job
        logger.info(f"Submitting slurm script {slurm_script} on {slurm_host}")
        command = (
            # f'MOOSE_WEBSERVERCONTROL_FILE_SOCKET="{socket_path}" '
            f"SLURMRUNNER_ARG=Controls/web_server/file_socket={socket_path} "
            f"sbatch {slurm_script}"
        )
        exit_code, result = self.exec_ssh_command(command)
        if exit_code != 0:
            self.error_exit(
                f"Failed to submit slurm script {slurm_script} "
                f"on {slurm_host}; result:\n\n{result}"
            )
        # Get the ID from submission on-screen output
        match_job = match(r"^Submitted batch job (\d+)$", result)
        if not match_job:
            raise RuntimeError(f"Failed to parse submitted slurm job from '{result}'")
        job_id = int(match_job.group(1))
        logger.info(f"Submitted job {job_id} on {slurm_host}")
        self._job_id = job_id

        # Wait for the job to start
        logger.info(f"Waiting for job {job_id} to start on {slurm_host}...")
        node = None
        while True:
            state, node = self.get_job_state()
            if state != "PENDING":
                if state == "COMPLETED":
                    self.error_exit(f"Job {job_id} on {slurm_host} has completed")
                elif state == "RUNNING":
                    logger.info(f"Job {job_id} on {slurm_host} running on {node}")
                    break
                else:
                    self.error_exit(
                        f"Unknown job state {state} for job {job_id} on {slurm_host}"
                    )

            sleep(1)

        # Wait for the socket to exist
        assert node is not None
        logger.info(f"Waiting for socket {node}:{socket_path} to exist...")
        while True:
            state, _ = self.get_job_state()
            if state != "RUNNING":
                self.error_exit(f"Job {job_id} on {slurm_host} has state {state}")
            exit_code, _ = self.exec_ssh_command(
                f"ssh -o StrictHostKeyChecking=no {node} test -e {socket_path}"
            )
            if exit_code == 0:
                logger.info(f"Socket {node}:{socket_path} exists")
                break

            sleep(1)

        # Setup the proxy to connect to the socket on the compute node
        logger.info(f"Setting up {node}:{socket_path} proxy through {self.slurm_host}")
        try:
            subprocess.run(
                [
                    "ssh",
                    "-q",
                    # Just forward, no remote shell
                    "-N",
                    # Background after authentication
                    "-f",
                    # Exit quickly if the forward fails
                    "-o",
                    "ExitOnForwardFailure=yes",
                    "-o",
                    "StreamLocalBindUnlink=yes",
                    # Disable SSH key checking
                    "-o",
                    "StrictHostKeyChecking=no",
                    "-o",
                    "ConnectTimeout=10",
                    "-o",
                    "ServerAliveInterval=1",
                    # Forward on remote to local
                    "-L",
                    f"{socket_path}:{socket_path}",
                    # Jump through the login node
                    "-J",
                    self.slurm_host,
                    # And onto the login node
                    node,
                ],
                # text=True,
                check=True,
                # stderr=subprocess.PIPE,
                # stdout=subprocess.PIPE,
            )
        except subprocess.CalledProcessError as e:
            print(e.stdout)
            raise ConnectionError(
                f"Failed to setup proxy to {node} through {self.slurm_host}"
            ) from e
        logger.info(f"Proxy to {node}:{socket_path} setup")

        SocketRunner.initialize(self, data)
