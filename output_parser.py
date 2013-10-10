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
		_rearrangment_breaker = False
		_junction_breaker = False
		_fields_breaker = False
		self.hits_v = [] # to be parsed in another function
		self.hits_d = [] # to be parsed in another function
		self.hits_j = [] # to be parsed in another function
		self.domain_classification = ""
		self.blast_dict = {}
		for line in query:
			if "Query" in line:
				self.query = line.split(":")[1].strip()
			if "Domain classification requested:" in line:
				self.domain_classification = line.split(":")[1].strip()
			if "rearrangement summary" in line:
				self.rearrangment_summary_titles = line.strip().split("(")[2].split(")")[0].split(",")
				_rearrangment_breaker = True
				continue
			if _rearrangment_breaker:
				self.rearrangment_summary = line.strip().split("\t")
				_rearrangment_breaker = False
			if "junction details" in line:
				self.junction_detail_titles = line.strip().split("(")[2].split(")")[0].split("\t")
				_junction_breaker = True
				continue
			if _junction_breaker:
				self.junction_detail = line.strip().split("\t")
				_junction_breaker = False
			if "Alignment summary" in line:
				self.alignment_summary_titles = line.strip().split("(")[1].split(")")[0].split(",")
				self.alignment_summary_titles = ["section"] + self.alignment_summary_titles
			if line.startswith("FR3-IMGT"):
				self.fr3_alignment_summary = line.strip().split()
			if line.startswith("CDR3-IMGT"):
				self.cdr3_alignment_sumary = [line.strip().split()[0]] + line.strip().split()[2:]
			if line.startswith("Total"):
				self.total_alignment_summary = line.strip().split()
			if "# Fields:" in line:
				self.hit_fields = line.strip().split(":")[1].split(",")
				_fields_breaker = True
			if _fields_breaker:
				if line.startswith("V"):
					self.hits_v.append(line)
				elif line.startswith("D"):
					self.hits_d.append(line)
				elif line.startswith("J"):
					self.hits_j.append(line)
		
		self.process()
		print self.blast_dict
		sys.exit()
		

	def process(self):
		self.blast_dict[self.query] = {"domain_classification":self.domain_classification,
										   "rearranment":self.parse_rearranment(),
										   "junction":self.parse_junction(),
										   "fr3_align":self.parse_fr3_align(),
										   "cdr3_align":self.parse_cdr3_align(),
										   "total_align":self.parse_total_align(),
										   "v_hits":self.parse_v_hits(),
										   "d_hits":self.parse_d_hits(),
										   "j_hits":self.parse_j_hits()
										   }
	def parse_rearranment(self):
		_return_dict = {}
		for title,value in zip(self.rearrangment_summary_titles,self.rearrangment_summary):
			if len(value.split(',')) > 1: 
				_return_dict[title] = tuple(value.split(',')) #cast multiple entries for tuple, makes them easier for json
			else:
				_return_dict[title] = value
		return _return_dict

	def parse_junction(self):
		_return_dict = {}
		for title,value in zip(self.junction_detail_titles,self.junction_detail):
			if "(" in value:
				_return_dict["d_or_j"] = value.split("(")[1].split(")")[0]
			else:
				_return_dict[title] = value
		return _return_dict

	def parse_fr3_align(self):
		_return_dict = {}
		for title,value in zip(self.alignment_summary_titles,self.fr3_alignment_summary):
			_return_dict[title] = value
		return _return_dict

	def parse_cdr3_align(self):
		_return_dict = {}
		for title,value in zip(self.alignment_summary_titles,self.cdr3_alignment_sumary):
			_return_dict[title] = value
		return _return_dict
	
	def parse_total_align(self):
		_return_dict = {}
		for title,value in zip(self.alignment_summary_titles,self.total_alignment_summary):
			_return_dict[title] = value
		return _return_dict

	def parse_v_hits(self):
		_return_dict = {}
		rank = 1
		for entry in self.hits_v:
			_return_dict['rank'] = rank
			for value,title in zip(entry.split()[1:],self.hit_fields):
				_return_dict[title] = value
			rank += 1
		return _return_dict

	def parse_d_hits(self):
		_return_dict = {}
		rank = 1
		for entry in self.hits_d:
			_return_dict['rank'] = rank
			for value,title in zip(entry.split()[1:],self.hit_fields):
				_return_dict[title] = value
		rank += 1
		return _return_dict

	def parse_j_hits(self):
		_return_dict = {}
		rank = 1
		for entry in self.hits_j:
			_return_dict['rank'] = rank
			for value, title in zip(entry.split()[1:],self.hit_fields):
				_return_dict[title] = value
		rank += 1
		return _return_dict


if __name__ == '__main__':
	igblast_output("igblast_out_all.blast_out")