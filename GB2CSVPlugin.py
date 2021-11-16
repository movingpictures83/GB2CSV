import sys

class GB2CSVPlugin:
    def input(self, outfile):
       self.gbfile = open(outfile, 'r')
    
    def run(self):
       self.header = ["locus_tag", "type", "+/-", "start", "end", "wstart", "wend", "inference", "note", "codon_start", "transl_table", "pseudo", "product", "protein_id", "translation", "gene", "EC_number", "anticodon", "ncRNA_class", "db_xref", "transl_except"]

       line=""
       while not line.startswith("FEATURES"):
          line = self.gbfile.readline()

       #print(line) # FEATURE line

       while not line.startswith("gene    "):
          line = self.gbfile.readline().strip()

       #print(line) # first gene line

       self.GBentries = dict()
       while (True):
        currentgene = ""
        while not (line.startswith("CDS") or line.startswith("rRNA") or line.startswith("tRNA") or line.startswith("ncRNA") or line.startswith("tmRNA")):
           line = self.gbfile.readline().strip()
           # Assumes at least one
           if (line.startswith("/locus_tag")):
               locus_tag = line[line.index('=')+1:].strip()
               currentgene = locus_tag
               self.GBentries[currentgene] = dict()

        #print("CURRENT GENE:"+currentgene)
        contents = line.split()
        #print(contents)
        entry = contents[1]
        #print("ENTRY: "+entry)
        self.GBentries[currentgene]['+/-'] = "+"
        if (line.startswith("CDS")):
            self.GBentries[currentgene]['type'] = "CDS"
        elif (line.startswith("rRNA")):
            self.GBentries[currentgene]['type'] = "rRNA"
        elif (line.startswith("tRNA")):
            self.GBentries[currentgene]['type'] = "tRNA"
        elif (line.startswith("ncRNA")):
            self.GBentries[currentgene]['type'] = "ncRNA"
        elif (line.startswith("tmRNA")):
            self.GBentries[currentgene]['type'] = "tmRNA"
        else:
            print("WARNING Unknown Type:"+line)
            exit()
        if (entry.startswith("complement")):
           entry = entry[len("complement")+1:len(entry)-1]  # Strip complement
           self.GBentries[currentgene]['+/-'] = "-"
        if (entry.startswith("join")):
           firstend = int(entry[entry.find('..')+2:entry.find(',')])
           secondend = int(entry[entry.rfind('..')+2:len(entry)-1])
           if (secondend > firstend):
              self.GBentries[currentgene]['start'] = entry[entry.find('(')+1:entry.find('.')]
              self.GBentries[currentgene]['end'] = str(secondend)
              self.GBentries[currentgene]['wstart'] = ""
              self.GBentries[currentgene]['wend'] = ""
           else:
              self.GBentries[currentgene]['start'] = entry[entry.find('(')+1:entry.find('.')]
              self.GBentries[currentgene]['end'] = entry[entry.find('..')+2:entry.find(',')]
              self.GBentries[currentgene]['wstart'] = entry[entry.find(',')+1:entry.rfind('..')]
              self.GBentries[currentgene]['wend'] = entry[entry.rfind('..')+2:len(entry)-1]
           pass
        else: # regular, no join
           self.GBentries[currentgene]['start'] = entry[entry.find('(')+1:entry.find('.')]
           self.GBentries[currentgene]['end'] = entry[entry.find('..')+2:]
           self.GBentries[currentgene]['wstart']=""
           self.GBentries[currentgene]['wend']=""

        lasttag = ""
        self.GBentries[currentgene]['inference'] = ""
        self.GBentries[currentgene]['note'] = ""
        self.GBentries[currentgene]['codon_start'] = ""
        self.GBentries[currentgene]['transl_table'] = ""
        self.GBentries[currentgene]['pseudo'] = "No"
        self.GBentries[currentgene]['product'] = ""
        self.GBentries[currentgene]['protein_id'] = ""
        self.GBentries[currentgene]['translation'] = ""
        self.GBentries[currentgene]['gene'] = ""
        self.GBentries[currentgene]['EC_number'] = ""
        self.GBentries[currentgene]['anticodon'] = ""
        self.GBentries[currentgene]['ncRNA_class'] = ""
        self.GBentries[currentgene]['db_xref'] = ""
        self.GBentries[currentgene]['transl_except'] = ""
        while (not (line.startswith("ORIGIN") or line.startswith("gene    "))):
           line = self.gbfile.readline().strip()
           if (line.startswith("ORIGIN")):
               print("DONE")
               break
           elif (line.startswith("/locus_tag")):
               pass
           elif (line.startswith("/inference")):
               self.GBentries[currentgene]['inference'] = line[line.index('=')+1:]
               lasttag = 'inference'
           elif (line.startswith("/note")):
               self.GBentries[currentgene]['note'] = line[line.index('=')+1:]
               lasttag = 'note'
           elif (line.startswith("/codon_start")):
               self.GBentries[currentgene]['codon_start'] = line[line.index('=')+1:]
               lasttag = 'codon_start'
           elif (line.startswith("/transl_table")):
               self.GBentries[currentgene]['transl_table'] = line[line.index('=')+1:]
               lasttag = 'transl_table'
           elif (line.startswith("/pseudo")):
               self.GBentries[currentgene]['pseudo'] = "Yes"
               lasttag = 'pseudo'
           elif (line.startswith("/product")):
               self.GBentries[currentgene]['product'] = line[line.index('=')+1:]
               lasttag = 'product'
           elif (line.startswith("/protein_id")):
               self.GBentries[currentgene]['protein_id'] = line[line.index('=')+1:]
               lasttag = 'protein_id'
           elif (line.startswith("/translation")):
               self.GBentries[currentgene]['translation'] = line[line.index('=')+1:]
               lasttag = 'translation'
           elif (line.startswith("/gene")):
               self.GBentries[currentgene]['gene'] = line[line.index('=')+1:]
               lasttag = 'gene'
           elif (line.startswith("/EC_number")):
               self.GBentries[currentgene]['EC_number'] = line[line.index('=')+1:]
               lasttag = 'EC_number'
           elif (line.startswith("/anticodon")):
               self.GBentries[currentgene]['anticodon'] = line[line.index('=')+1:]
               lasttag = 'anticodon'
           elif (line.startswith("/ncRNA_class")):
               self.GBentries[currentgene]['ncRNA_class'] = line[line.index('=')+1:]
               lasttag = 'ncRNA_class'
           elif (line.startswith("/db_xref")):
               self.GBentries[currentgene]['db_xref'] = line[line.index('=')+1:]
               lasttag = 'db_xref'
           elif (line.startswith("/transl_except")):
               self.GBentries[currentgene]['transl_except'] = line[line.index('=')+1:]
               lasttag = 'transl_except'
           elif (not line.startswith("gene    ")):
               self.GBentries[currentgene][lasttag] += " " + line
           if (line.startswith("gene    ")):
               pass
               #print(line)
               #print("Going back now")
           elif (lasttag != "" and self.GBentries[currentgene][lasttag].count(',') != 0):
               self.GBentries[currentgene][lasttag] = self.GBentries[currentgene][lasttag].replace(',', ' ')
        if (line.startswith("ORIGIN")):
           break
 #print(self.GBentries[currentgene])


#for line in self.gbfile:
#    if (line.startswith("gene")):
#        while (not 

    def output(self, filename):
       outfile = open(filename, 'w')
       for i in range(len(self.header)):
           outfile.write(self.header[i])
           if (i == len(self.header)-1):
               outfile.write("\n")
           else:
               outfile.write(',')

       for key in self.GBentries:
           outfile.write(key+",")
           print(key)
           print(self.GBentries[key])
           for i in range(1, len(self.header)): # Don't do locus_tag
               outfile.write(self.GBentries[key][self.header[i]])
               if (i == len(self.header)-1):
                   outfile.write("\n")
               else:
                   outfile.write(',')


