import requests
import os

class GetRest(object):

    def __init__(self, function):
        self.url_base = ''
        self.function = function
        self.headers = None

    def get_header(self):
        self.headers = {'Content-Type': 'application/json',
           'Authorization': 'Bearer ' + os.environ['token_api']}
    
    def get_url_base(self):
        self.url_base = os.environ['endpoint_api']

    def get_new_authorization(self, error):
        if error == 401:
            print('Need to create a new token')
        elif error == 403:
            print('Login data is not vaild')
        else:
            print('Error number %d', error)

    def performRequest(self):
        self.get_header()
        self.get_url_base()
        url_get = self.url_base + self.function
        result_request = requests.get(url_get, headers=self.headers)
        if result_request.status_code == 200:
            return result_request.json()
        else:
            self.get_new_authorization(result_request.status_code)