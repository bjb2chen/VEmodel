import re
phoneRegex = re.compile(r'\(\d\d\d\) \d\d\d-\d\d\d\d')
resume = '''Benny Jun Chen
bjb2chen@uwaterloo.ca | (647) 382-1693
EDUCATION
UNIVERSITY OF WATERLOO
CANDIDATE FOR
MASTER OF SCIENCE,
THEORETICAL CHEMISTRY
January 2022 - Present
BACHELOR OF SCIENCE,
HONOURS CHEMISTRY,
COMPUTATIONAL SPECIALIZATION,
CO-OPERATIVE PROGRAM
September 2016 - April 2021
Waterloo, Ontario, Canada
COURSEWORK
Quantum Chemistry
Computational Chemistry
Inorganic Chemistry
Organic Chemistry
Analytical Chemistry
Python Programming
Computer Hardware and Systems
Business Ethics and HR
Philosophy of Critical Thinking
TECHNICAL SKILLS
Python
Jupyter Notebook
Gaussian
QuantumESPRESSO
SAP and ERP Software
Inventory and SDS Databases
UV Materials
SOCIETIES AND CLUBS
2019 | UW ChemClub | Member
2018 | Tespa | Competitor
2017 | UW League Club | Member
2017 | UW Hearthstone Club | Member
2016 | UW Science Society | Member
PERSONAL INTERESTS
Beverage Science • Geopolitics • Reddit •
eSports • Swimming • Personal Finance •
Cats • Cryptography • String Theory •
Minimalism • Settlers of Catan
LANGUAGES
English (native) • Cantonese (native) •
Mandarin (beginner fluency)
RELEVANT WORK EXPERIENCE
NATURAL RESOURCES CANADA | CHEM ENG/CHEMISTRY
January 2020 - August 2020 | Ottawa, ON
• Developed a dynamic model in Python to simulate the separation performance
ofrare earth elements via electrodialysis.
• Administered bench-scale electrodialysis experiments to evaluate the efficacy
of separating rare earth elements.
• Supported CanmetMINING’s research efforts in the Rare Earth Elements and
Chromite RD Program.
UNIVERSITY OF WATERLOO | CHEMICAL COMPLIANCE ASSISTANT
September 2018 - December 2018 | Waterloo, ON
• Managed University of Waterloo’s ’eRPortal’ 100,000+ chemical inventory
database forthe duration of the Fall 2018 term on behalf of UW Safety Office.
• Optimized the transition to WHMIS 2015.
• Introduced the solution of using internal campus Xerox printers forlabel
printing instead of external label printers.
• Pursued outreach with more than 75 labs campus-wide, amiably providing
training, and completed inventories totaling 7,000 chemicals.
AXALTA COATING SYSTEMS | ASSOCIATE CHEMIST
January 2018 - April 2018 | Cornwall, ON
• Formulated industrial wood coatings in the RnD laboratory by applying
understanding of chemistry and dispersion dynamics.
• Specialized in UV solvents,resins, powders, and photo-initiators: unilaterally
produced 75 percent of all UV formulations in the lab.
• Performed chemical tests and analyses to assist company salesmen concerns.
• Collaborated with internal manufacturing processes between quality control,
shipping, and technical production plant.
CANADA COMPUTERS | HEAD OFFICE RECEPTIONIST |
CUSTOMER SERVICE LEAD REPRESENTATIVE
June 2016 - August 2016, May 2017 - August 2017 | Richmond Hill, ON
• Obtained advance promotion to Receptionist in recognition of exceeding
manager expectations and stellar performance relative to peers, lead of 4
member phone support team.
• Increased clients served by 17 percent through phone and online inquiry form.
• Resolved customerinquiries competently via developing ad-hoc compromising
solutions that satisfied both clients and management.
ADDITIONAL EXPERIENCE
SCINAPSE UNDERGRADUATE COMPETITION | EXECUTIVE LEADER
October 2017 | Waterloo, ON
• Exhibited leadership by overseeing a team of 3 other aspiring scientists
towards researching practical solutions to hydraulic fracturing and formulation
of a successful case study.
'''
mo = phoneRegex.findall(resume)
print(mo)

digitRegex = re.compile(r'(0|1|2|3|4|5|6|7|8|9)') # is equivalent to jut r'\d', shorthand char class

lyrics = '''12 drummers drumming,
11 pipers piping,
10 lords a-leaping,
9 ladies dancing,
8 maids a-milking,
7 swans a-swimming,
6 geese a-laying,
5 gold rings,
4 calling birds,
3 French hens,
8 turtle doves,
And 1 partridge in a pear tree!
'''
xmasRegex = re.compile(r'\d+ \w+')
print(xmasRegex.findall(lyrics))

vowelRegex = re.compile(r'[aeiouAEIOU]') # r'(a|e|i|o|u)''
print(vowelRegex.findall('Robocop eats baby food.'))
doubleVowelRegex = re.compile(r'[aeiouAEIOU]{2}') 
print(doubleVowelRegex.findall('Robocop eats baby food.'))

consonantsRegex = re.compile(r'[^aeiouAEIOU]') # caret makes it the negative case
print(consonantsRegex.findall('Robocop eats baby food.'))