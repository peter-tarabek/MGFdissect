-- this script updates the values that exist in the second table ("ions" in this example)
UPDATE ions
SET
	blank = (SELECT blank 
             FROM res1
             WHERE ions.ID = res1.compound_id);