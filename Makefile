profile:
	python3 -m cProfile -s 'time' clipkit-runner.py -i test_files/12_YIL115C_Anc_2.253_aa_aln.fasta -o output/test > output/profile.txt

run:
	python3 -m clipkit-runner -i test_files/12_YIL115C_Anc_2.253_aa_aln.fasta -o output/test

install:
	# install so clipkit command is available in terminal
	python setup.py install

develop:
	# https://setuptools.readthedocs.io/en/latest/setuptools.html#development-mode
	python setup.py develop

test:
	python -m pytest

test_default_params:
	# test clipkit with default parameters
	python3 -m clipkit-runner -i test_files/test.fa -o output/test.fa.test_default_params
	python3 -m clipkit-runner -i test_files/12_YIL115C_Anc_2.253_codon_aln.fasta -o output/12_YIL115C_Anc_2.253_codon_aln.fasta.test_default_params
	python3 -m clipkit-runner -i test_files/12_YIL115C_Anc_2.253_aa_aln.fasta -o output/12_YIL115C_Anc_2.253_aa_aln.fasta.test_default_params
	cmp --silent test_expected_output/test.fa.test_default_params output/test.fa.test_default_params || echo "Fail"
	cmp --silent test_expected_output/12_YIL115C_Anc_2.253_codon_aln.fasta.test_default_params output/12_YIL115C_Anc_2.253_codon_aln.fasta.test_default_params || echo "Fail"
	cmp --silent test_expected_output/12_YIL115C_Anc_2.253_aa_aln.fasta.test_default_params output/12_YIL115C_Anc_2.253_aa_aln.fasta.test_default_params || echo "Fail"
