# Piecewise-regression
This script requires input in the form of constants and an input text file.
Constants include:
1. Latitude and longitude of study site in digital form (default = Santa Rita Experimental Range, Arizona).
2. Number of hives (number of columns - 1).
3. Number of iterations for the segmented() function (default = 100).
4. Number of samples per hour (default = 12, or every 5 minutes).
5. Output file name.

Input file needs to have the following format:
1. First column is date and time of sample, with the format MM/DD/YYYY hh:mm
2. All other columns are weight values. Column headers will appear in output file.

Included is a sample input file.
