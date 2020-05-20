#!/usr/bin/env python
import psycopg2
connection = psycopg2.connect(user = "meta_user",
                          password = "",
                          host = "localhost",
                          port = "5432",
                          database = "string_sample_tables")

curs = connection.cursor()
conn = psycopg2.connect(user = "meta_user",
                              password = "",
                              host = "localhost",
                              port = "5432",
                              database = "string_tables")
cursor = conn.cursor()

curs.execute("SELECT string_id FROM refseq_to_string")
fetch = curs.fetchall()
for group in fetch:
	cursor.execute("SELECT * FROM protein_interaction_data WHERE protein_a='%s'"%group)
	fetched = cursor.fetchall()
	if fetched != None:
		for groups in fetched:
			curs.execute("INSERT INTO protein_interaction_data(protein_a,protein_b,mode,action_is_directional,action_is_acting,score,extra) VALUES(%s,%s,%s,%s,%s,%s,%s)",groups)
curs.close()
cursor.close()
connection.commit()
conn.commit()
