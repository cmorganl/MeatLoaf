# MeatLoaf

In this directory are a collection of scripts I have developed to perform relatively simple weightlifting. The majority are small scripts that I've judged to be used more than once but by very few people.

There is no organization - think of it like one of those unkept thrift shops in the cool part of town. Maybe its all a load of junk... but maybe you'll find something that will work.

Welcome to MeatLoaf. It has a little bit of everything.

# Script usage

## JGI\_portal\_downloader.py

This script can be used for programmatically downloading all project files (or a subset thereof) from the JGI's download portal using their API. It just requires a user’s credentials and the portal name for the sample(s). Details on how this works can be found at https://genome.jgi.doe.gov/portal/help/download.jsf#/api.

To create a 'cookies' file, log in using the following command:
```
curl 'https://signon.jgi.doe.gov/signon/create' --data-urlencode 'login=USER_NAME' --data-urlencode 'password=USER_PASSWORD' -c cookies > /dev/null
```

The steps for collecting these portal-specific sample names are as follows:
1. Go to https://genome.jgi.doe.gov/portal/
2. Search for (meta)genome name
3. Click “Download” for each row under the “Resources” column
4. Copy the name between the second and third ‘/’, usually the string after “portal/”.

To get the usage and options, run
```
python JGI_portal_downloader.py -h
```

