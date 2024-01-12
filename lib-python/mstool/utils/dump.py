import sqlite3
import pandas as pd

def dumpsql(ifile):
	prefix = '.'.join(ifile.split('.')[:-1])
	fmt    = ifile.split('.')[-1]
	assert fmt.lower() == 'dms', 'DMS format'
	conn = sqlite3.connect(ifile)

	with open(prefix + '.sql', 'w') as f:
		for line in conn.iterdump():
			f.write('%s\n' %line)

def dumpdf(ifile, table):
	prefix = '.'.join(ifile.split('.')[:-1])
	fmt    = ifile.split('.')[-1]
	assert fmt.lower() == 'dms', 'DMS format'
	conn = sqlite3.connect(ifile, 
		isolation_level=None, 
		detect_types=sqlite3.PARSE_COLNAMES)

	df = pd.read_sql_query("SELECT * FROM {:s}".format(table), conn)
	return df



