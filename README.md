## File descriptions
- `table_S1.xlsx`, `table_S3.xlsx`, and `table1.xlsx` are transcribed from the corresponding tables in the original paper
- `cohn_reader.R` reads the Excel files, changes data formats, summarizes data, and checks for consistency
- `replication.R` has code for replicating the original study. It's dependent on output from `cohn_reader.R`
- `theoretical_simul.R` has code for simulating the theoretical processes related to the data, to show why the original analyses were inappropriate
