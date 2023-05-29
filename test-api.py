from clipkit import clipkit
from clipkit.files import FileFormat

trim_run = clipkit(
    # raw_alignment=">1\nA-GTAT\n>2\nA-G-AT\n>3\nA-G-TA\n>4\nAGA-TA\n>5\nACa-T-\n",
    input_file_path="tests/integration/samples/simple.fa",
    # output_file_path="./programmatic-temp.phylip",
    # output_file_format=FileFormat.phylip,
    # mode=TrimmingMode.smart_gap,
    mode="smart-gap",
    gaps=None,
    sequence_type="nt",
)

print(trim_run.site_classification_counts)
