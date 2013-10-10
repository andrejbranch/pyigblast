import sys

class igblast_output():
	'''Pass the hole output to this class and it will deal with it
	Args: 
	file - string for the filename to parse
	output - string for the filename to the output
	'''
	def __init__(self,file):
		breaker = True
		query_holder = []
		for line in open(file):
			if "IGBLASTN" in line:
				breaker = False
				if query_holder:
					single_blast_entry(query_holder)
				query_holder = []
				continue
			if not breaker:
				query_holder.append(line)	
		single_blast_entry(query_holder)

class single_blast_entry():
	'''The helper class to parse an individual blast result'''
	def __init__(self,query):
		for line in query:
			if "Query" in line:
				self.query = line.split(":")[1].strip()
				print self.query
if __name__ == '__main__':
	igblast_output("igblast_out_all.blast_out")