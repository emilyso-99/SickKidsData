#!/usr/bin/env python

import psycopg2
import sys
import decimal
connection = psycopg2.connect(user = "meta_user",
                              password = "",
                              host = "localhost",
                              port = "5432",
                              database = "string_tables")

curs = connection.cursor()
# READING FROM INPUT FILE 
# f = open('/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/input_data.txt','r')
# curs.copy_from(f,'input_data')
# f.close()
# READING FROM INPUT FILE 

#making refseq to string
# curs.execute("SELECT refseq_id,ncbi_taxid FROM input_data")
# result = curs.fetchall()
# # this for loop collects the string id 
# with open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/input_data.txt","w") as outfile:
#   outfile.write("Refseq ID" + "\t" + "NCBI tax ID" + "\n")
#    # + "String ID" + "\t" + "Cluster ID" + "\n")
#   for refseq_id in result:
#     refseq = (refseq_id[0]).split(".")
#     outfile.write(refseq[0] + "\t" + refseq_id[1] + "\n")
# with open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/input_data.txt","r") as infile, open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/sample_refseq_to-st.txt","w") as outfile:
#   outfile.write("refseq_id" + "\t" + "ncbi_taxid\tstring_id\n")
#   outfile.flush()
#   next(infile)
#   for line in infile:
#     line_split = line.split("\t")
#     ncbi_taxid = (line_split[1].split("\n"))[0]
#     refseq_id = line_split[0]
#     curs.execute("SELECT string_id FROM refseq_to_string WHERE refseq_id='%s'" % refseq_id)
#     fetched = curs.fetchone()
#     if (fetched != None):
#       print("Got Through")
#       sys.stdout.flush()
#       outfile.write(refseq_id + '\t' + ncbi_taxid + '\t' + fetched[0] + '\n')
#       outfile.flush()
# with open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/input_data.txt","r") as infile, open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/finding_clusters_info.txt","w") as outfile:
#   outfile.write("ncbi_taxid" + "\t" + "cluster_id\tdescription" + "\n")
#   outfile.flush()
#   next(infile)
#   for line in infile:
#     line_split = line.split('\t')
#     ncbi_taxid = (line_split[1].split("\n"))[0]
#     curs.execute("SELECT cluster_id,description FROM stringdb_clusters_info WHERE ncbi_taxid='%s'" % ncbi_taxid)
#     fetched = curs.fetchall()
#     for fetch in fetched:
#       if (fetch != None):
#         print("Got Through")
#         sys.stdout.flush()
#         outfile.write(ncbi_taxid + '\t' + fetch[0] + '\t' + fetch[1] + '\n')
#         outfile.flush()
# with open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/sample_refseq_to_string.txt","r") as infile, open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/string_gene_ontology","w") as outfile:
#   outfile.write("string_id\t" + "ncbi_taxid\t" + "gene_ontology_id" + "\t" + "gene_ontology_category" + "\n")
#   outfile.flush()
#   next(infile)
#   for line in infile:
#     line_split = line.split("\t")
#     string_id = (line_split[1].split("\n"))[0]
#     ncbi_taxid = line_split[0]
#     curs.execute("SELECT * FROM stringdb_gene_ontology WHERE ncbi_taxid='%s'" % ncbi_taxid)
#     fetched = curs.fetchall()
#     if (fetched != None):
#       for fetch in fetched:
#         if (string_id == fetch[3]):
#           outfile.write(string_id+'\t'+ ncbi_taxid + '\t' + fetch[2]+'\t'+fetch[1] + '\n')
#           outfile.flush()
conn = psycopg2.connect(user = "meta_user",
                              password = "",
                              host = "localhost",
                              port = "5432",
                              database = "eggnog_tables")

cursor = conn.cursor()
# with open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/refseq_to_string.txt","r") as infile, open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/finding_cluster_proteins.txt","w") as outfile:
#   outfile.write("string_protein" + "\t" + "cluster_id" + "\t" + "ncbi_taxid\n")
#   outfile.flush()
#   next(infile)
#   for line in infile:
#     line_split = line.split("\t")
#     string_id = (line_split[2].split("\n"))[0]
#     ncbi_taxid = line_split[1]
#     curs.execute("SELECT cluster_id FROM stringdb_clusters_proteins WHERE protein_member='%s'" % string_id)
#     fetched = curs.fetchall()
#     for fetch in fetched:
#       if (fetch != None):
#         print("Got Through")
#         sys.stdout.flush()
#         outfile.write(string_id + '\t' + fetch[0] + '\t' + ncbi_taxid + '\n')
#         outfile.flush()
# with open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/string_to_enog.txt","r") as infile, open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/finding_bactnog_proteins.txt","w") as outfile:
#   outfile.write("string_protein" + ":" + "enog_id" + ":cog_category:description:" + "proteinds\n")
#   outfile.flush()
#   next(infile)
#   for line in infile:
#     line_split = line.split("\t")
#     description = (line_split[3].split("\n"))[0]
#     group_id = line_split[1]
#     string_protein = line_split[0]
#     cogcategory = line_split[2]
#     cursor.execute("SELECT proteinids FROM eggnog_bactnog_members WHERE cluster_id='%s'" % group_id)
#     fetched = cursor.fetchone()
#     if (fetched != None):
#       print("Got Through")
#       sys.stdout.flush()
#       outfile.write(string_protein + ':' + group_id + ':' + cogcategory + ':' + description + ':' + fetched[0] + '\n')
#       outfile.flush()
# with open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/sample_refseq_to_string.txt","r") as file_1, open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/proteins_gene_ontology.txt","w") as file_2:
#   file_2.write("ncbi_taxid\tstringid\tgo_id\tgocategory\n")
#   file_2.flush()
#   next(file_1)
#   for line in file_1:
#     line_split = line.split("\t")
#     string_id = (line_split[2].split("\n"))[0]
#     curs.execute("SELECT gene_ontology_id,go_category,ncbi_taxid FROM stringdb_gene_ontology WHERE string_id='%s'" % string_id)
#     fetch = curs.fetchone();
#     if fetch != None:
#       ncbi_taxid = fetch[2]
#       file_2.write(ncbi_taxid + '\t' + string_id + '\t' + fetch[0]+'\t'+fetch[1]+'\n')
#       file_2.flush()
# with open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/sample_refseq_to_string.txt","r") as infile_1, open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/refseq_to_enog.txt","w") as outfile:
#   outfile.write("string_id\tgroupid\tcog_category\tdescription\n")
#   outfile.flush()
#   next(infile_1)
#   for line in infile_1:
#     line_split = line.split('\t')
#     refseq_id = line_split[0]
#     string_id = (line_split[2].split("\n"))[0]
#     curs.execute("SELECT uniprot FROM uniprot_refseq WHERE refseq='%s'" % refseq_id)
#     uniprot = curs.fetchone()
#     if uniprot != None:
#       cursor.execute("SELECT luca_og,bac_og FROM uniprotenog WHERE uniprot_id='%s'" % uniprot[0])
#       ogs = cursor.fetchone()
#       if ogs != None:
#         print("OG:%s" % ogs[1])
#         if (ogs[1] == "-"):
#           print(string_id,ogs[0])
#           print("Nope")
#           sys.stdout.flush()
#           outfile.write(string_id + '\t' + ogs[0] + '\tN/A\tN/A\n')
#           outfile.flush()
#         else:
#           sys.stdout.flush()
#           cursor.execute("SELECT cog_category,description FROM eggnog_bactnog_annotations WHERE group_id='%s'" % ogs[1])
#           cog = cursor.fetchone()
#           if cog != None:
#             print
#             print(cog[0],cog[1])
#             outfile.write(string_id + '\t' + ogs[1] + '\t' + cog[0] + '\t' + cog[1] + '\n')
#             outfile.flush()
# with open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/string_to_enog.txt","r") as infile_1, \
# open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/string_to_enog_with_interactions.txt","w") as outfile:
#   outfile.write("ncbi taxid\trefseq id\tstring_id\tgroupid\tcog_category\tinteracting protein\tinteracting protein cog \n")
#   outfile.flush()
#   next(infile_1)
#   for line in infile_1:
#     line_split = line.split('\t')
#     string_id = line_split[0]
#     group_id = line_split[1]
#     cog_category = line_split[2]
#     description = line_split[3]
#     curs.execute("SELECT DISTINCT protein_b FROM stringdb_protein_interactions WHERE protein_a='%s'" % string_id)
#     fetched = curs.fetchall()
#     curs.execute("SELECT ncbi_taxid,refseq_id FROM refseq_to_string WHERE string_id='%s'"%string_id)
#     refseq_tax = curs.fetchone()
#     if fetched != None:
#       for fetch in fetched:
#           call = "SELECT cluster_id FROM eggnog_bactnog_members WHERE proteinids LIKE" +"'%" + fetch[0] + "%'"
#           cursor.execute(call)
#           cluster = cursor.fetchone()
#           if cluster != None:
#             outfile.write(refseq_tax[0] + '\t' + refseq_tax[1] + '\t'+ string_id + '\t' + group_id + '\t' + cog_category + '\t' + fetch[0] +  '\t' + cluster[0] + '\n')
#           else:
#             outfile.write(refseq_tax[0] + '\t' + refseq_tax[1] + '\t'+ string_id + '\t' + group_id + '\t'+ cog_category + '\t' + fetch[0] +  '\t' + "N/A" + '\n')
#     else:
#       outfile.write(refseq_tax[0] + '\t' + refseq_tax[1] + '\t'+ string_id+'\t'+group_id+'\t'+cog_category+'\t'+"N/A\t"+"N/A\n")
# with open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/string_to_enog_with_interactions.txt","r") as file_1, \
# open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/interaction_output_data.txt","w'") as output:
#   output.write("refseq_id\tstring_id\tgroupid\tcog_category\tinteracting protein\tinteracting protein cog \n")
#   output.flush()
#   next(file_1)
#   for line in file_1:
#     line_split = line.split("\t")
#     string_id = line_split[0]
#     curs.execute("SELECT refseq_id FROM refseq_to_string WHERE string_id='%s'" % string_id)
#     fetch = curs.fetchone();
#     if fetch != None:
#       refseq_id = fetch[0]
#       output.write(refseq_id+'\t'+line)
#       output.flush()
# with open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/string_to_enog.txt","r") as file_1, \
# open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/string_to_enog_with_go.txt","w'") as output:
#   output.write("string_id\tgroupid\tcog_category\tdescription\tgo_id\tgo_category\n")
#   output.flush()
#   next(file_1)
#   for line in file_1:
#     line_split = line.split("\t")
#     cluster_id = line_split[1]
#     string_id = line_split[0]
#     cursor.execute("SELECT proteinids FROM eggnog_bactnog_members WHERE cluster_id='%s'" % cluster_id)
#     fetch = cursor.fetchone();
#     if fetch != None:
#       fetch_split = (fetch[0]).split(",")
#       string = None
#       for proteinid in fetch:
#         curs.execute("SELECT go_category,gene_ontology_id FROM stringdb_gene_ontology WHERE string_id='%s'"%proteinid)
#         go_terms = curs.fetchone()
#         print(go_terms)
#         if go_terms != None:
#           string = "Output"
#           output.write((line.split('\n'))[0] + '\t' + go_terms[1] + '\t' + go_terms[0] + '\n')
#           output.flush()
#       print(string)
#       if (string == None):
#         output.write(line)
#         output.flush()
# with open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/string_to_enog.txt","r") as file_1, \
# open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/taxonomic_ratio.txt","w'") as output:
#   output.write("groupid\tcog_category\tdescription\tgo_id\tgo_category\n")
#   output.flush()
#   next(file_1)
#   for line in file_1:
#     line_split = line.split("\t")
#     cluster_id = line_split[1]
#     string_id = line_split[0]
#     cursor.execute("SELECT proteinids FROM eggnog_bactnog_members WHERE cluster_id='%s'" % cluster_id)
#     fetch = cursor.fetchone();
#     if fetch != None:
#       fetch_split = (fetch[0]).split(",")
#       for proteinid in fetch:
#         curs.execute("SELECT go_category,gene_ontology_id FROM stringdb_gene_ontology WHERE string_id='%s'"%proteinid)
#         go_terms = curs.fetchone()
#         print(go_terms)
#         if go_terms != None:
#           output.write((line.split('\n'))[0] + '\t' + go_terms[1] + '\t' + go_terms[0] + '\n')
#           output.flush()
#           break
#       else:
#         output.write(line)
#         output.flush

# with open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/string_to_enog_with_interactions.txt","r") as file,\
# open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/counting.txt","w") as out:
# 	out.write("cluster\tcount\n")
# 	for line in file:
# 		cluster = (line.split('\t'))[3]
# 		cursor.execute("SELECT proteincount FROM eggnog_bactnog_members WHERE cluster_id='%s'" % cluster)
# 		cluster_id = cursor.fetchone()
# 		if cluster_id != None:
# 			out.write(cluster+'\t'+str(cluster_id[0])+'\n')
# with open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/Emily_Stable2.txt","r") as file_2, \
# open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/new_input_data.txt","w") as out:
# 	curs.execute("CREATE TABLE new_input_data(refseq varchar,ncbi_taxid varchar,locus varchar,name varchar,strain varchar,genus varchar,supergroup varchar)")
# 	for line in file_2:
# 		line_split = line.split('\t')
# 		description = line_split[1].split('|')
# 		if len(description) >= 9:
# 			locus = description[8]
# 			name = description[6]
# 		else:
# 			locus = "N/A"
# 			name = description[4]
# 		strain = line_split[4]
# 		genus = line_split[5]
# 		supergroup = line_split[9]
# 		refseq = description[3]
# 		ncbi_taxid = line_split[3]
# 		out.write(refseq+'\t'+ncbi_taxid+'\t'+locus+'\t'+name+'\t'+strain+'\t'+genus+'\t'+supergroup+'\n')
# 		out.flush()
# with open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/new_input_data.txt","r") as out:
# 	curs.copy_from(out,'new_input_data')
# 	out.close()
# with open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/counting.txt","r") as file,\
# open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/edited_my_taxonomy.txt","r") as file_2,\
# open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/family_tax_ratios.txt","w") as out:
# 	out.write("cluster\tcount\n")
# 	out.flush()
# 	species_dict = {}
# 	for line in file:
# 		cluster = (line.split('\t'))[0]
# 		cursor.execute("SELECT proteinids FROM eggnog_bactnog_members WHERE cluster_id='%s'" % cluster)
# 		cluster_id = cursor.fetchone()
# 		if cluster_id != None:
# 			proteinids = cluster_id[0].split(',')
# 			for protein in proteinids:
# 				curs.execute("SELECT ncbi_taxid,refseq_id FROM refseq_to_string WHERE string_id='%s'" % protein)
# 				ids = curs.fetchone()
# 				if ids != None:
# 					ncbi_taxid = ids[0]
# 					cursor.execute("SELECT taxid_lineage FROM eggnog_taxonomic_classification WHERE ncbi_taxid='%s'" % ncbi_taxid)
# 					lineage = cursor.fetchone()
# 					if lineage != None:
# 						for line in file_2:
# 								if (line.split('\t'))[0] in lineage:
# 							if (line.split('\t'))[1] == "family":
# 									if cluster not in species_dict:
# 										species_dict[cluster] = []
# 									species_dict[cluster].append(protein)
# 	for key in species_dict:
# 		out.write(key+'\t'+species_dict[key]+'\n')
# 		out.flush()
# with open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/string_to_enog.txt","r") as inputfile,\
# open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/final_output_info.txt",'w') as outputfile:
# 	outputfile.write("locus\tname\trefseq\tncbi_taxid\tstrain\tgenus\tsupergroup\n")
# 	for line in inputfile:
# 		line_split = line.split('\t')
# 		string_id = line_split[0]
# 		curs.execute("SELECT refseq_id FROM refseq_to_string WHERE string_id='%s'"% string_id)
# 		fetch = curs.fetchone()
# 		if fetch != None:
# 			sql_statement = "SELECT * FROM new_input_data WHERE refseq LIKE '%" + fetch[0] + "%'"
# 			curs.execute(sql_statement)
# 			info = curs.fetchone()
# 			if info != None:
# 				outputfile.write(info[2]+'\t'+info[3]+'\t'+info[0]+'\t'+info[1]+'\t'+info[4]+'\t'+info[5]+'\t'+info[6]+'\t' + line)
# 				outputfile.flush()
# with open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/finding_bactnog_proteins.txt","r") as countingfile,\
# open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/members_by_tax_id.txt",'w') as outputfile:
#       outputfile.write("cluster\ttaxidmembers\n")
#       outputfile.flush()
#       next(countingfile)
#       for line in countingfile:
#             family_string = ''
#             line_split = line.split(':')
#             cluster = line_split[1]
#             cursor.execute("SELECT proteincount FROM eggnog_bactnog_members WHERE cluster_id='%s'" % cluster)
#             count = cursor.fetchone()
#             if count != None:
#                   proteincount = str(count[0])
#                   print("Protein_count:"+proteincount+'\n')
#                   sys.stdout.flush()
#             proteins = ((line_split[4]).split('\n'))[0]
#             protein_split = proteins.split(',')
#             for protein in protein_split:
#                   curs.execute("SELECT ncbi_taxid FROM refseq_to_string WHERE string_id='%s'"%protein)
#                   tax_id = curs.fetchone()
#                   if tax_id != None:
#                         family_string += tax_id[0] + ","
#                         print(len(family_string.split(',')))
#                         sys.stdout.flush()
#             outputfile.write(cluster+'\t'+family_string+'\n')
#             outputfile.flush()
# with open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/finding_bactnog_proteins.txt","r") as countingfile,\
# open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/members_by_tax_id.txt",'w') as outputfile:
#       outputfile.write("cluster\ttaxidmembers\n")
#       outputfile.flush()
#       next(countingfile)
#       for line in countingfile:
#             family_string = ''
#             line_split = line.split(':')
#             cluster = line_split[1]
#             cursor.execute("SELECT proteincount FROM eggnog_bactnog_members WHERE cluster_id='%s'" % cluster)
#             count = cursor.fetchone()
#             if count != None:
#                   proteincount = str(count[0])
#                   print("Protein_count:"+proteincount+'\n')
#                   sys.stdout.flush()
#             proteins = ((line_split[4]).split('\n'))[0]
#             protein_split = proteins.split(',')
#             for protein in protein_split:
#                   curs.execute("SELECT ncbi_taxid FROM refseq_to_string WHERE string_id='%s'"%protein)
#                   tax_id = curs.fetchone()
#                   if tax_id != None:
#                         family_string += tax_id[0] + ","
#                         print(len(family_string.split(',')))
#                         sys.stdout.flush()
#             outputfile.write(cluster+'\t'+family_string+'\n')
#             outputfile.flush()

# with open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/family_ratios.txt","w") as outfile,\
# open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/new_members_by_tax_id.txt","r") as infile:
# 	outfile.write("eggnog_id\ttaxid\tnumber\n")
# 	outfile.flush()
# 	genus_dict = {}
# 	next(infile)
# 	for line in infile:
# 		line_split = line.split('\t')
# 		cluster = line_split[0]
# 		if cluster not in genus_dict:
# 			genus_dict[cluster] = {}
# 		proteinids = ((line_split[1]).split('\n'))[0]
# 		proteins = proteinids.split(',')
# 		print(cluster)
# 		sys.stdout.flush()
# 		for protein in proteins:
# 			cursor.execute("SELECT taxid_lineage FROM eggnog_taxonomic_classification WHERE ncbi_taxid='%s'"%protein)
# 			lineage = cursor.fetchone()
# 			if lineage != None:
# 				protein_list = (lineage[0]).split(',')
# 				print("PROTEIN:%s"%protein)
# 				cursor.execute("SELECT ncbi_taxid FROM direct_taxonomic_levels WHERE rank='family'")
# 				fetched = cursor.fetchall()
# 				if fetched != None:
# 					for fetch in fetched:
# 						if fetch[0] in protein_list:
# 							print("FETCH:%s"%fetch)
# 							sys.stdout.flush()
# 							if fetch[0] not in genus_dict[cluster]:
# 								genus_dict[cluster][fetch[0]] = 0
# 							genus_dict[cluster][fetch[0]] += 1
# 		print("now for dictionary")
# 		sys.stdout.flush()
# 		for key in genus_dict[cluster]:
# 			print("CLUSTER%s"%cluster)
# 			sys.stdout.flush()
# 			print("KEY:%s"%key)
# 			sys.stdout.flush()
# 			outfile.write(cluster+'\t'+key+'\t'+str(genus_dict[cluster][key])+'\n')
# 			outfile.flush()
with open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/genus_ratios.txt","r") as infile,\
open("/Users/emilyso/Documents/ThirdYear/BCB330/microbial_metafunctional_maps/EmilySo/final_genus_ratios.txt","w") as outfile:
	outfile.write("eggnog_id\ttaxid\tratio\n")
	outfile.flush()
	genus_dict = {}
	next(infile)
	for line in infile:
		line_split = line.split('\t')
		cluster = line_split[0]
		tax = line_split[1]
		count = (line_split[2]).split('\n')[0]
		cursor.execute("SELECT proteincount FROM eggnog_bactnog_annotations WHERE group_id='%s'"%cluster)
		counts = cursor.fetchone()
		if counts!= None:
			print("CLUSTER:%s"%cluster)
			print("COUNT:%s"%str(count[0]))
			sys.stdout.flush()
			outfile.write(cluster+'\t'+tax+'\t'+str(decimal.Decimal(count)/decimal.Decimal(counts[0]))+'\n')
			outfile.flush()
cursor.close()
conn.commit()
conn.close()
connection.commit()
curs.close()
connection.close()






