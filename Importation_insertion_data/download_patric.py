import re
from robobrowser import RoboBrowser

# Browse to Rap Genius
browser = RoboBrowser(history=True)
browser.open('https://www.patricbrc.org')


# Search for Queen
form = browser.get_form(action='/dijit__WidgetsInTemplateMixin_6/')

form["username"] = 'gresch'
form["password"] = 'sequencing'
browser.session.headers['Referer'] = base_url

browser.submit_form(form)
print(str(browser.select))