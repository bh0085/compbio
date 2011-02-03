import xml 

# 3 handler functions

def start_element(name, attrs):
    print 'Start element:', name, attrs
def end_element(name):
    print 'End element:', name
def char_data(data):
    print 'Character data:', repr(data)

p = xml.parsers.expat.ParserCreate()

p.StartElementHandler = start_element
p.EndElementHandler = end_element
p.CharacterDataHandler = char_data
p.Parse(open('ARCH.xml').read(), 1)


from xml.dom.minidom import parse, parseString

dom1 = parse('ARCH.xml') # parse an XML file by name


raise Exception()
