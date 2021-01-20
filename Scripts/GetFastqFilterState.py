#!/usr/bin/env python

import pandas as pd
import json
import re
import os
import sys

## defined fetech information
def getSampleState(injson, sample=None):
    with open(injson, 'r') as js:
        data = json.load(js)
    del data['summary']['before_filtering']['read1_mean_length']
    del data['summary']['before_filtering']['read2_mean_length']
    del data['summary']['after_filtering']['read1_mean_length']
    del data['summary']['after_filtering']['read2_mean_length']
    del data['filtering_result']['too_long_reads']
    del data['filtering_result']['passed_filter_reads']

    dd1 = pd.DataFrame([data['summary']['before_filtering']])
    dd1.columns = ["_".join(["before",i]) for i in dd1]

    dd2 = pd.DataFrame([data['filtering_result']])

    dd3 = pd.DataFrame([data['summary']['after_filtering']])
    dd3.columns = ["_".join(["after",i]) for i in dd3]

    ddd = dd1.join(dd2)
    ddd = ddd.join(dd3)
    if sample:
        ddd.index = [sample]

    return(ddd)

## output
if __name__ == "__name__":
    RE_SAMPLE = re.compile(r'Step01\.FastqFilter\/(\S+)\/(\S+)\.json')
    out = pd.DataFrame()

    with open(sys.argv[1], 'r') as jsf:
        for jsfile in jsf.readlines():
            jsfile = jsfile.rstrip()
            sample = RE_SAMPLE.search(str(jsfile)).group(1)
            abc = getSampleState(jsfile, sample)
            out = out.append([abc])
        res = pd.DataFrame(out)
    
    res.to_csv(sys.argv[2], sep="\t")