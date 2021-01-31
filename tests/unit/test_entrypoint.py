import os
import pytest
import subprocess


class TestEntrypoint(object):
    def test_help(self):
        cmd = "clipkit --help"
        exit_status = os.system(cmd)
        assert exit_status == 0

    def test_run(self):
        cmd = "clipkit tests/integration/samples/simple.fa"
        exit_status = os.system(cmd)
        assert exit_status == 0

    def test_input_error(self):
        cmd = "clipkit /file/doesnt/exist"
        response = subprocess.check_output([cmd], stderr=subprocess.STDOUT, shell=True)
        assert response == b"Input file does not exist\n"

    def test_run_no_args(self):
        cmd = "clipkit"
        exit_status = os.system(cmd)
        assert exit_status == 0
