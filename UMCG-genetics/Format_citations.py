from pybtex.database.input import bibtex

#open a bibtex file
parser = bibtex.Parser()
bibdata = parser.parse_file("citations.bib")

#loop through the individual references
Miss = []
records = []
for bib_id in bibdata.entries:
    b = bibdata.entries[bib_id].fields
    try:
        # change these lines to create a SQL insert
        Title = b["title"]
        if not 'journal' in b:
            Jorunal = 'Preprint'
        else:
            Jorunal = b["journal"]
        Year = b["year"]
        if not 'volume' in b:
            Volume = ""
        else: 
            Volume = b['volume']
        if not 'number' in b:
            number = "" 
        else:
            number = b['number']
        if not 'pages' in b:
            pages = ""
        else: 
            pages = b['pages']
        #deal with multiple authors
        N_author = len(bibdata.entries[bib_id].persons["author"])
        Authors = []
        if N_author <= 10:
            for author in bibdata.entries[bib_id].persons["author"]:
                Name = ''.join(author.first()) + ' ' + ''.join(author.last())
                Name = Name.replace(r"{\'a}", "a").replace(r"{\'A}", "A").replace(r"{\'o}", "o").replace(r"{\`e}", "e").replace(r"{\"o}", "o")
                if ''.join(author.first()) == "Sergio": Name = Name + "*"
                Authors.append(Name)
        else:
            for author in bibdata.entries[bib_id].persons["author"][0:7]:
                Name = ''.join(author.first()) + ' ' + ''.join(author.last())
                Name = Name.replace(r"{\'a}", "a").replace(r"{\'A}", "A").replace(r"{\'o}", "o").replace(r"{\`e}", "e").replace(r"{\"o}", "o")
                if ''.join(author.first()) == "Sergio": Name = Name + "*"
                Authors.append(Name)
            Authors.append('et al. ')
            #for author in bibdata.entries[bib_id].persons["author"][-3:N_author]:
            #    Name = ''.join(author.first()) + ' ' + ''.join(author.last())
            #    Name = Name.replace('{\'a}', "a")
            #    Name = Name.replace('{\'A}', "A")
            #    Authors.append(Name)
        Fields ='. '.join([",".join(Authors), Year, Title, Jorunal, f'{Volume}.{number}', pages])
        records.append({"Year": int(Year), "Fields": Fields})
    # field may not exist for a reference
    except KeyError as e:
        Miss.append({"bib_id": bib_id, "error": str(e)})
        continue
sorted_records = sorted(records, key=lambda x: x["Year"], reverse=True)
for record in sorted_records:
    with open('formatted_citations.txt', 'a' ) as F:
        F.write(record["Fields"]+"\n")

# Print missing records with errors
if Miss:
    print("\nMissing records:")
    for missing in Miss:
        print(f"Bib ID: {missing['bib_id']} | Error: {missing['error']}")