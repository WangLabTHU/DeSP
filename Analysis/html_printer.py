'''
Date: 2020-07-06 14:11:44
LastEditTime: 2021-06-23 19:42:20
LastEditors: Lekang
Organization: BBNC Lab@THU
Describtions: Do not edit
'''
def color_bold(s,color = '#ff1111'):
    return f'<b style = "color: {color}">{s}</b>'

def background_color_bold(s,color = '#ff1111'):
    return f'<b style = "background-color: {color}; color: black;">{s}</b>'

def html_templete(templete_name,path = 'Analysis/templete.html'):
    with open(path) as f:
        out = ''
        # locate desired paragraph
        while True:
            line = f.readline()
            if len(line) >= 4 and line[:2]== '<!':
                if line[4:-4] == templete_name:
                    break
            if line == '':
                return None
        
        # read the paragraph
        while True:
            line = f.readline()
            if line == '' or line[:2] == '<!': break
            else: out += line

        return out

class html_table:
    def __init__(self,templete = 'table_style'):
        self.style = html_templete(templete)

    def print(self,table,title,styled = True):
        out = ''
        # print title
        out += self.print_row(title, head = True)
        for row in table:
            out += self.print_row(row)
        out = '<table>\n' + out + '</table>\n'
        if styled: out = self.style + out
        return out

    def print_row(self,row,head = False):
        out = ''
        
        if type(row) is tuple:
            assert len(row) == 2
            tr_style, row = row
        else:
            tr_style, row = '', row
        
        for element in row:
            if type(element) is tuple:
                assert len(element) ==  2
                element, style = element
            else:
                element, style = element, ''
            
            if head:
                out += '\t\t<th ' + style + '>' + str(element) + '</th>\n'
            else:
                out += '\t\t<td ' + style + '>' + str(element) + '</td>\n'
           
        out = '\t<tr ' + tr_style + '>\n' + out + '\t</tr>\n'  
        return out