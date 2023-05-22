from clipkit.api import clipkit
from clipkit.modes import TrimmingMode

result = clipkit(
    raw_alignment=">1\nA-GTAT\n>2\nA-G-AT\n>3\nA-G-TA\n>4\nAGA-TA\n>5\nACa-T-\n",
    # input_file_path="tests/integration/samples/simple.fa",
    # output_file_path="./programmatic-temp.fa",
    output_file_format=None,
    mode=TrimmingMode.smart_gap,
    gaps=None,
)

print(result)
