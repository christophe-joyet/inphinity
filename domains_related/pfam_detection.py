# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 15:27:26 2017

@author: Stage
"""

from __future__ import print_function
from six.moves import urllib

from files_treatment.domain_json_file_hmm import *

from objects_new.temp_prot import *

import logging

import json

class SmartRedirectHandler(urllib.request.HTTPRedirectHandler):
    def http_error_302(self, req, fp, code, msg, headers):
        logging.debug(headers)
        return headers



class PFAM_detection(object):
    def __init__(self, sequence, id_prot = -1):
        self.sequence = '>No name \n' + sequence
        self.base_url = 'https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan'
        self.json_file = None
        self.id_prot = id_prot

        
    def detecterPFAM(self):
        sequence = self.sequence
        handle = self.hmmscan(hmmdb = 'pfam', seq = sequence)
        self.json_file = json.loads(handle.read().decode('utf8'))
        
        json_domain = DomainJsonFile(self.json_file)
        json_domain.search_domains()
        #insert id_protein
        temp_prot = id_prot_temp()
        temp_prot.id_prot_temp = self.id_prot

        temp_prot.create_id_prot()


        return json_domain.list_domains
        

    def _hmmer(self, endpoint, args1, args2):
        opener = urllib.request.build_opener(SmartRedirectHandler())
        urllib.request.install_opener(opener);
    
        params = urllib.parse.urlencode(args1).encode("utf-8")
        try:
            req = urllib.request.Request(endpoint,
                                  data = params,
                                  headers={"Accept" : "application/json"})
            v = urllib.request.urlopen(req, timeout = 60)
        except urllib.error.HTTPError as e:
            raise Exception("HTTP Error 400: %s" % e.read())
        except timeout:
            print("Timeout error")
    
        results_url = v['location']
    
        enc_res_params = urllib.parse.urlencode(args2)
        modified_res_url = results_url + '?' + enc_res_params
    
        results_request = urllib.request.Request(modified_res_url)
        f = urllib.request.urlopen(results_request)
        return f

    def phmmer(self, **kwargs):
        """Search a protein sequence against a HMMER sequence database.
        Arguments:
          seq - The sequence to search -- a Fasta string.
          seqdb -- Sequence database to search against.
          range -- A string range of results to return (ie. 1,10 for the first ten)
          output -- The output format (defaults to JSON).
        """
        logging.debug(self, kwargs)
        args = {'seq' : kwargs.get('seq'),
                'seqdb' : kwargs.get('hmmdb')}
        args2 = {'output' : kwargs.get('output', 'json'),
                 'range' : kwargs.get('range', None)}
        return self._hmmer("https://www.ebi.ac.uk/Tools/hmmer/search/phmmer", args, args2)
    
    def hmmscan(self, **kwargs):
        logging.debug(kwargs)
        args = {'seq' : kwargs.get('seq'),
                'hmmdb' : kwargs.get('hmmdb')}
        args2 = {'output' : 'json'}
        range = kwargs.get('range', None)
        if range:
            args2['range'] = range
        return self._hmmer("https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan", args, args2)
