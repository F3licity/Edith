# import urllib
# url='http://bs-solutions.nl/temp/crypt.html'
# urllib.unquote(url).decode('utf8')
# u'example.com?title=\u043f\u0440\u0430\u0432\u043e\u0432\u0430\u044f+\u0437\u0430\u0449\u0438\u0442\u0430'
# print urllib.unquote(url).decode('utf8')
# example.com?title=hallo+hallo

str = 'aHR0cDovL2JzLXNvbHV0aW9ucy5ubC90ZW1wL2NyeXB0Lmh0bWw='
# str.decode(encoding='UTF-8',errors='strict')

# str = "this is string example....wow!!!";
str = str.encode('UTF-8','strict');

print "Encoded String: " + str;
print "Decoded String: " + str.decode('UTF-8','strict')