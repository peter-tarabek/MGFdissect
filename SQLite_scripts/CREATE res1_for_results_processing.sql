-- this script creates a mirror of the base table with columns for found compounds 
-- in multiple samples (in columns)
CREATE TABLE res1 AS
SELECT compound_id, "HighestIon-m/z", "RT/min", ESIpolarity, ionType, Compound, CAS, Formula, EN_use -- adjust the selection of columns as needed
FROM base;
-- now add the columns...
ALTER TABLE res1 ADD COLUMN "sample1" INTEGER;
ALTER TABLE res1 ADD COLUMN "sample2" INTEGER;
ALTER TABLE res1 ADD COLUMN "sample3" INTEGER;
ALTER TABLE res1 ADD COLUMN "sample4" INTEGER;
ALTER TABLE res1 ADD COLUMN "sample5" INTEGER;
ALTER TABLE res1 ADD COLUMN "sample6" INTEGER;
ALTER TABLE res1 ADD COLUMN "sample7" INTEGER;
ALTER TABLE res1 ADD COLUMN "sample8" INTEGER; -- more sample names (columns) can be added
ALTER TABLE res1 ADD COLUMN blank INTEGER;