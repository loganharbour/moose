# This file is part of the MOOSE framework
# https://mooseframework.inl.gov
#
# All rights reserved, see COPYRIGHT for full restrictions
# https://github.com/idaholab/moose/blob/master/COPYRIGHT
#
# Licensed under LGPL 2.1, please see LICENSE for details
# https://www.gnu.org/licenses/lgpl-2.1.html

"""Defines the SlurmRunner."""

import paramiko
import os

from logging import getLogger
from typing import Optional, Tuple
from re import match

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

    def exec_ssh_command(self, command: str, check: bool = True) -> Tuple[int, str]:
        client = self.get_ssh_client()
        try:
            _, stdout, stderr = client.exec_command(command)
        except Exception as e:
            raise RuntimeError(
                f"Failed to run SSH command {command} on {self.slurm_host}"
            ) from e
        exit_code = stdout.channel.recv_exit_status()
        result = "".join(stdout.readlines())
        if exit_code != 0:
            result += "".join(stderr.readlines())
        return exit_code, result.rstrip()

    def initialize(self, data: dict):
        slurm_host = self.slurm_host
        slurm_script = self.slurm_script

        # Setup the SSH client
        self.get_ssh_client()

        # Make sure the slurm path exists
        if self.exec_ssh_command(f"test -f {slurm_script}")[0] != 0:
            raise FileNotFoundError(
                f"Slurm script {slurm_script} on {slurm_host} does not exist"
            )

        # Submit the job
        logger.info(f"Submitting slurm script {slurm_script} on {slurm_host}")
        command = f"sbatch {slurm_script}"
        exit_code, result = self.exec_ssh_command(command)
        if exit_code != 0:
            print(result)
            raise RuntimeError(
                f"Failed to submit slurm script {slurm_script} on {slurm_host}"
            )
        match_job = match(r"^Submitted batch job (\d+)$", result)
        if not match_job:
            raise RuntimeError(f"Failed to parse submitted slurm job from '{result}'")
        job_id = int(match_job.group(1))
        logger.info(f"Submitted job {job_id} on {slurm_host}")

        raise SystemExit("foo")

        SocketRunner.initialize(self, data)
