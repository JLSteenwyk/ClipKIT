profile:
	# python3 -m cProfile -s 'time' clipkit-runner.py test_files/12_YIL115C_Anc_2.253_aa_aln.fasta -o output/test > output/profile.txt
	python3 -m cProfile -s 'time' clipkit-runner.py  tests/integration/samples/EOG092C0CZK_aa_aln.fasta -o output/test > output/large_profile.txt

run:
	python3 -m clipkit-runner test_files/12_YIL115C_Anc_2.253_aa_aln.fasta -o output/test --log

run.simple:
	python3 -m clipkit-runner test_files/test.fa -o output/test --log	

install:
	# install so clipkit command is available in terminal
	python setup.py install

develop:
	# https://setuptools.readthedocs.io/en/latest/setuptools.html#development-mode
	python setup.py develop

test: test.unit test.integration

test.unit:
	python -m pytest -m "not integration"

test.integration:
	python -m pytest --basetemp=output -m "integration"

test_default_params:
	# TODO: finish migrating to tests/integration/test_gappy_mode.py
	# TODO: Jacob -- write more integration tests for different functions, etc.
	# test clipkit with default parameters
	python3 -m clipkit-runner test_files/test.fa -o output/test.fa.test_default_params
	python3 -m clipkit-runner test_files/12_YIL115C_Anc_2.253_codon_aln.fasta -o output/12_YIL115C_Anc_2.253_codon_aln.fasta.test_default_params
	python3 -m clipkit-runner test_files/12_YIL115C_Anc_2.253_aa_aln.fasta -o output/12_YIL115C_Anc_2.253_aa_aln.fasta.test_default_params
	cmp --silent test_expected_output/test.fa.test_default_params output/test.fa.test_default_params || echo "Fail"
	cmp --silent test_expected_output/12_YIL115C_Anc_2.253_codon_aln.fasta.test_default_params output/12_YIL115C_Anc_2.253_codon_aln.fasta.test_default_params || echo "Fail"
	cmp --silent test_expected_output/12_YIL115C_Anc_2.253_aa_aln.fasta.test_default_params output/12_YIL115C_Anc_2.253_aa_aln.fasta.test_default_params || echo "Fail"
