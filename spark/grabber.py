import requests, re, random, time

# year, month, day_start, day_end
tasks = [
    ['2017','05','0100','3120'],
    ['2017','04','0100','3020'],
    ['2017','03','0100','3120'],
    ['2017','02','0100','2820'],
    ['2017','01','0100','3120'],
    ['2016','12','0100','3120'],
    ['2016','11','0100','3020'],
    ['2016','10','0100','3120'],
    ['2016','09','0100','3020'],
    ['2016','08','0100','3120'],
    ['2016','07','0100','3120'],
    ['2016','06','0100','3020'],
    ['2016','05','0100','3120'],
    ['2016','04','0100','3020'],
    ['2016','03','0100','3120'],
    ['2016','02','0100','2920'],
    ['2016','01','0100','3120'],
    ['2015','12','0100','3120'],
    ['2015','11','0100','3020'],
    ['2015','10','0100','3120'],
    ['2015','09','0100','3020'],
    ['2015','08','0100','3120'],
    ['2015','07','0100','3120'],
    ['2015','06','0100','3020']
]

time_pattern = re.compile('(\d+)\/(\d+)')

all_outlines = []

for task in tasks:
    print(task)
    
    r = requests.get('http://www.weather.uwyo.edu/cgi-bin/sounding?region=europe&TYPE=TEXT%3ALIST&YEAR='+task[0]+'&MONTH='+task[1]+'&FROM='+task[2]+'&TO='+task[3]+'&STNM=12425')
    lines = [line for line in r.content.decode('utf-8').splitlines()]
    del lines[:4]

    output = False
    outlines = []
    measure = []
    for i, line in enumerate(lines):
        if i >= 6 and '<H2>' in lines[i-6]:
            output = True
        elif 'Station information' in line:
            output = False
        elif 'Observation time' in line:
            datetime = time_pattern.findall(line)[0]
            if '0000' not in datetime[1]:
                dt_col = datetime[0] + datetime[1]
                outlines.extend([dt_col+'\t'+row for row in measure])
                measure = []
        if output:
            cols = [col.strip() for col in line.split(' ') if col.strip() != '']
            if len(cols) == 11:
                measure.append('\t'.join(cols[:5]))
        
    all_outlines.extend(outlines)
    
fs = open('data.txt', 'wb')
fs.write('\n'.join(all_outlines).encode('utf-8'))
fs.close()