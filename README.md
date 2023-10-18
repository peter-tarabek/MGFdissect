# MGFdissect  
## Short description  
*"dissecting" and reordering MGF*
   
The tool in this branch parses and groups MGF (mascot generic format) MS/MS spectra from multiple files into a single file. The reordering and grouping of particular entries (between BEGIN IONS and END IONS statements) is based on their parent ion m/z and retention time values. A minimalistic user interface is provided through svDialogs package.     
This script apart from combining the MFGs searches also a custom .db database (SQLite, the search table must be named "ions", example included in the repo as "nts_db.db") and assigns a numeric value based on previous analysis relating to the presence of compounds in blanks (0_ identified previously in samples but not present in blank, 1_ identified previously and also present in blank, NA_ present in the database but not identified previously).  
It also reduces the number of MGF items (spectra) - saves only the most intense ones in the final merged file.  
