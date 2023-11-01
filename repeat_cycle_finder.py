import subprocess

# Loop the number of times you want to run memote
for i in range(0, 10):
    # Define the output file name
    output_file = "report_" + str(i) + ".html"

    # Define the command as a list of strings
    command = [
        "memote",
        "report",
        "snapshot",
        "--filename",
        output_file,
        "model.xml",
    ]

    # Run the command using subprocess
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Wait for the process to finish
    stdout, stderr = process.communicate()
