import copy
class weightinfo:
    id=''
    pdf=''
    muF=''
    muR=''
    name=''

    def __init__(self):
        self.id=''
        self.pdf=''
        self.muF=''
        self.muR=''
        self.name=''
                                                                                                                        
class parser:
    
    _filename=''#=open("v242_LO_weight_info.txt",'r')
    
    _startgroup='' #='<weightgroup'
    _endgroup='' #'</weightgroup>'
    _startweight=''
    _endweight=''
    ket='>'   
    bra='<'
    Group_list=[]
    mode=0
    title=''
    def __init__(self):
        self._filename=''#=open("v242_LO_weight_info.txt",'r')
        self._startgroup='' #='<weightgroup'
        self._endgroup='' #'</weightgroup>'
        self._startweight=''
        self._endweight=''
        self.ket='>'
        self.bra='<'
        self.Group_list=[]
        self.mode=0
        self.title=''
                                        
        
    def define_braket(self,bra,ket):
        self.ket=ket
        self.bra=bra
    class Group:
        
        title=''
        isscale=0
        weights=[]
        #idlist=[]
        #pdflist=[]
        #def add_id(self,idx):
        #    self.idlist.append(idx)
        #def add_pdf(self,pdf):
        #    self.pdflist.append(pdf)
        def __init__(self):
            self.title=''
            self.weights=[]
        def set_title(self,title):
            self.title=title
        def size(self):
            return len(self.weights)
        def add_weight(self, weight):
            self.weights.append(weight)
        def clear(self):
            del self.weights[:]
            #self.title=''

                                                
    def set_filename(self,name):
        self._filename=name
    def set_startgroup(self,name):
        self._startgroup=name
    def set_endgroup(self,name):
        self._endgroup=name
    def set_startweight(self,name):
        self._startweight=name
    def set_endweight(self,name):
        self._endweight=name
                                
    def title_parser(self,title):
        #print "title="+title
        sptitle=title.split()
        #print sptitle
        mode=0
        newtitle=''
        for i in range(len(sptitle)):
            if "name=" in sptitle[i]:
                mode=1
            if mode==1:
                newtitle=newtitle+sptitle[i]

                #        if newtitle.endswith(self.ket):
                #            newtitle=newtitle[:-1]
                #if newtitle.startswith('name="') and newtitle.endswith('"'):
        if newtitle.startswith('name="'):
            newtitle=newtitle.lstrip('"name=')
            #newtitle=newtitle.rstrip('"')
            newtitle=newtitle[0:newtitle.find('"')]
            #            print newtitle
        return newtitle




    def weight_parser(self,weight):
        #print weight
        weight=weight.replace("PDF=    ","PDF=") ##For v242NLO
        weight=weight.replace("PDF=   ","PDF=") ##For v242NLO
        weight=weight.replace("PDF=  ","PDF=") ##For v242NLO
        weight=weight.replace("PDF= ","PDF=") ##For v242NLO
        #print "weight="+weight
        spweight=weight.split()
        #        print spweight
        #print len(spweight)
        id=''
        pdf=''
        muF=''
        muR=''
        name=''
                                                
        
        mode=0 ##mode1=>weight value mode2=>weight name
        pdfmode=0
        #print spweight
        for part in spweight:
            
            # print part+str(mode)
            if self._startweight in part: mode=1
            if self._endweight in part:
                name+=part
                name=name.rstrip(self._endweight)
                mode=0

            if mode==1:

                if ("MUF=" in part) or ("muF" in part):
                    #print "muf"
                    muF=part
                    muF=muF.lstrip('muF=')
                    muF=muF.lstrip('MUF="')
                    muF=muF.rstrip('"')
                    muF=muF.rstrip(self.ket)
                if ("MUR=" in part) or ("muR" in part):
                    muR=part
                    muR=muR.lstrip('muR=')
                    muR=muR.lstrip('MUR="')
                    muR=muR.rstrip('"')
                    muR=muR.rstrip(self.ket)
                if "PDF=" in part and pdfmode==0:
                    #print "PDF matched="+part
                    pdf=part
                    pdf=pdf.lstrip('PDF="')
                    pdf=pdf.rstrip('"')
                    pdf=pdf.rstrip(self.ket)
                    pdfmode=1
                if "id=" in part:
                    id=part
                    id=id.lstrip('id="')
                    id=id.rstrip('"')
                    id=id.rstrip(self.ket)
                    id=id.rstrip('"')
#                if mode==1 and (self.ket in part):
#                    mode=2
#            elif mode==2:
#                print part+"!!!is after ket"
#                name+=part
        name=name.lstrip(self.ket)
        
        new_weight=weightinfo()
        new_weight.muF=muF
        new_weight.muR=muR
        new_weight.pdf=pdf
 
        new_weight.id=id
        new_weight.name=name
        return new_weight
        
    def run(self):
        f=open(self._filename,'r')
        #print "==="+self._filename+"======"
        lines=f.readlines()
        mode=0

        for line in lines:            
            #            print line
            #print "line="+line
            if mode==3 and ( self._startgroup in line ) :
                
                self.Group_list.append(tempGroup)
                mode=0
                
            if self._startgroup in line:
                
                
                title=line
                newtitle=self.title_parser(title)
                tempGroup=self.Group()
                
                tempGroup.set_title(newtitle)
                if "scale" in tempGroup.title:
                    tempGroup.isscale=1
                    
                mode=1

            elif self._endgroup in line:
#                newGroup=copy.deepcopy(tempGroup)
#                newGroup.weights=copy.deepcopy(tempGroup.weights)
#                self.Group_list.append(newGroup)
                self.Group_list.append(tempGroup)
                #print self.Group_list[0].title
#                print newGroup.title
 #               print len(newGroup.weights)

#                print str(self.Group_list[0].size())
#                tempGroup.clear()
                
                #print "after clear="+self.Group_list[0].title
                #print "after clear="+str(self.Group_list[0].size())
                #print "after clear="+newGroup.title
                #print "after clear="+str(len(newGroup.weights))
                mode=0

            elif mode==1:
                _weightinfo_=self.weight_parser(line)
                if _weightinfo_.id !='':tempGroup.add_weight(_weightinfo_)
                
            if mode==0 and (self._startweight in line):
                newtitle="unknown"
                tempGroup=self.Group()
                tempGroup.set_title(newtitle)
                mode=3 #without <weightgroup case
            if mode==3:
                _weightinfo_=self.weight_parser(line)
                if _weightinfo_.id !='':tempGroup.add_weight(_weightinfo_)
                #print  _weightinfo_.id

            if mode==3 and lines[-1]==line  : ##if last line
                
                self.Group_list.append(tempGroup)
                mode=0
                                            
#        for group in self.Group_list:
 #           print group.title
def runLO242():
    LO242=parser()
    LO242.set_filename('v242_LO_weight_info.txt')
    LO242.set_startgroup('<weightgroup')
    LO242.set_endgroup('</weightgroup>')
    LO242.set_startweight('<weight')
    LO242.set_endweight('</weight>')
    LO242.run()
    Group_list=LO242.Group_list
    # print len(Group_list)
    Ntotal=0
    #    for i in range(1,1080):
    #        mode=0
    #        for group in Group_list:
    #            print "======"+group.title+"=========Nweight="+str(len(group.weights))+"======"
    #            Ntotal+=len(group.weights)
    #            for weight in group.weights:
    #            print 'name="'+weight.name+'" id='+weight.id+" pdf="+weight.pdf+" muF="+weight.muF+" muR="+weight.muR
    #            print weight.id
    #                if i == int(weight.id)   : mode=1
    #        if mode==0 : print i
    #    print "Ntotal="+str(Ntotal)
    return LO242
def runLO26x():
    LO26x=parser()
    #Group_list2=LO26x.Group_list
    #print len(Group_list2)
    LO26x.define_braket('&lt;','&gt;')
    LO26x.set_filename('v26x_LO_weight_info.txt')
    LO26x.set_startgroup('&lt;weightgroup')
    LO26x.set_endgroup('&lt;/weightgroup&gt;')
    LO26x.set_startweight('&lt;weight')
    LO26x.set_endweight('&lt;/weight&gt;')
    LO26x.run()
    Group_list2=LO26x.Group_list
    #print len(Group_list2)
    Ntotal=0
    #    for i in range(1080,2116):
    #   mode=0
#    for group2 in Group_list2:
#        print "======"+group2.title+"=========Nweight="+str(len(group2.weights))+"======"
        #Ntotal+=len(group2.weights)
#        for weight2 in group2.weights:
                #print 'name="'+weight2.name+'" id='+weight2.id+" pdf="+weight2.pdf+" muF="+weight2.muF+" muR="+weight2.muR
#                print '" id='+weight2.id+" pdf="+weight2.pdf+" muF="+weight2.muF+" muR="+weight2.muR
#                print "Ntotal="+str(Ntotal)
#            if int(weight2.id)==i : mode=1
#        if mode==0 : print i
    return LO26x
def runNLO242():
    NLO242=parser()
    NLO242.set_filename('242NLO.log')
    NLO242.set_startgroup('<weightgroup')
    NLO242.set_endgroup('</weightgroup>')
    NLO242.set_startweight('<weight')
    NLO242.set_endweight('</weight>')
                    
    NLO242.run()
    Group_list=NLO242.Group_list

    return NLO242
def runNLO26x():
    NLO26x=parser()
    NLO26x.define_braket('&lt;','&gt;')
    NLO26x.set_filename('26xNLO.log')
    NLO26x.set_startgroup('&lt;weightgroup')
    NLO26x.set_endgroup('&lt;/weightgroup&gt;')
    NLO26x.set_startweight('&lt;weight')
    NLO26x.set_endweight('&lt;/weight&gt;')                    
    NLO26x.run()
    Group_list=NLO26x.Group_list
    
    return NLO26x
                                    
def compare_LO():
    LO242=parser()
    LO242.set_filename('v242_LO_weight_info.txt')
    LO242.set_startgroup('<weightgroup')
    LO242.set_endgroup('</weightgroup>')
    LO242.set_startweight('<weight')
    LO242.set_endweight('</weight>')
    LO242.run()
    Group_list=LO242.Group_list

    LO26x=parser()
    LO26x.define_braket('&lt;','&gt;')
    LO26x.set_filename('v26x_LO_weight_info.txt')
    LO26x.set_startgroup('&lt;weightgroup')
    LO26x.set_endgroup('&lt;/weightgroup&gt;')
    LO26x.set_startweight('&lt;weight')
    LO26x.set_endweight('&lt;/weight&gt;')
    LO26x.run()
    Group_list2=LO26x.Group_list

    '''    
    print "##Find matched weights"

    for group in Group_list:
        mode=0
        #print "======"+group.title+"=========Nweight="+str(len(group.weights))+"======"
        for group2 in Group_list2:
         #   print "======"+group2.title+"=========Nweight="+str(len(group2.weights))+"======"
            if group.title in group2.title or group2.title in group.title:mode=1
        if mode==1: print group.title+" Nweight="+str(len(group.weights)) 

    print "#####Find unmatched weights#######"
    for group in Group_list:
        mode=0
        #print "======"+group.title+"=========Nweight="+str(len(group.weights))+"======"
        for group2 in Group_list2:
            #   print "======"+group2.title+"=========Nweight="+str(len(group2.weights))+"======"
            if group.title in group2.title or group2.title in group.title:mode=1
        if mode==0: print group.title+" Nweight="+str(len(group.weights))

    for group2 in Group_list2:
        mode=0
        #print "======"+group.title+"=========Nweight="+str(len(group.weights))+"======"
        for group in Group_list:
            #   print "======"+group2.title+"=========Nweight="+str(len(group2.weights))+"======"
            if group.title in group2.title or group2.title in group.title:mode=1
        if mode==0: print group2.title+" Nweight="+str(len(group2.weights))
   '''         
    ####Check weight  v242 but not in v26x
    print "Check weight  v242 but not in v26x"
    N242=0
    for group in Group_list:
        print group.title+" Nweight="+str(len(group.weights))
        for weight1 in group.weights:
            N242+=1
            mode=0
            for group2 in Group_list2:
                for weight2 in group2.weights:
                    if int(weight1.pdf) == int(weight2.pdf) and float(weight1.muF)==float(weight2.muF) and float(weight1.muR)==float(weight2.muR) :
                        mode=1
                        #print weight1.id

            if mode==0 : print weight1.id
    print "N242="+str(N242)


    
    print "Check weight  v26x but not in v242"
    N26x=0
    for group2 in Group_list2:
        print group2.title+" Nweight="+str(len(group2.weights))
        for weight2 in group2.weights:
            N26x+=1
            mode=0
            for group in Group_list:
                for weight1 in group.weights:
                    if int(weight1.pdf) == int(weight2.pdf) and float(weight1.muF)==float(weight2.muF) and float(weight1.muR)==float(weight2.muR) :
                        mode=1
                        #print weight1.id
                        
            if mode==0 : print weight2.id
    print "N26x="+str(N26x)
    
    print "N26x-N242="+str(N26x-N242)
def List242():
    LO242=parser()
    LO242.set_filename('v242_LO_weight_info.txt')
    LO242.set_startgroup('<weightgroup')
    LO242.set_endgroup('</weightgroup>')
    LO242.set_startweight('<weight')
    LO242.set_endweight('</weight>')
    LO242.run()
    Group_list=LO242.Group_list
    for group in Group_list:
        #        print group.title+" pdf1="+group.weights[0].pdf+" #pdf="+str(len(group.weights))
#        print " pdf1="+group.weights[0].pdf+" #pdf="+str(len(group.weights))
        print group.weights[0].pdf+","

def List26x():
    LO26x=parser()
    LO26x.define_braket('&lt;','&gt;')
    LO26x.set_filename('v26x_LO_weight_info.txt')
    LO26x.set_startgroup('&lt;weightgroup')
    LO26x.set_endgroup('&lt;/weightgroup&gt;')
    LO26x.set_startweight('&lt;weight')
    LO26x.set_endweight('&lt;/weight&gt;')
    LO26x.run()
    Group_list2=LO26x.Group_list
    for group2 in Group_list2:
        print group2.title+" pdf1="+group2.weights[0].pdf+" #pdf="+str(len(group2.weights))
class group_template:
    start_pdf=-1
    end_pdf=-1
    npdf=0
    isscale=0
    title=''
    def __init__(self):
        self.start_pdf=-1
        self.npdf=0
        self.isscale=0
    def set_start_pdf(self,pdf):
        self.start_pdf=pdf
    def set_npdf(self,npdf):
        self.npdf=npdf
    def set_isscale(self,isscale):
        self.isscale=isscale
    def set_title(self,title):
        self.title=title
    def set_end_pdf(self,pdf):
        self.end_pdf=pdf
                    
def define_set(v242_group_list, centralpdf=''):##make template
    groups=[]

    for v242group in v242_group_list:
        group=group_template()
        #pdf_list.append(group.weights[0].pdf)
        group.set_start_pdf(v242group.weights[0].pdf)
        if v242group.weights[0].pdf=='' and bool(v242group.isscale) :
             group.set_start_pdf(centralpdf)
        group.set_npdf(len(v242group.weights))
        group.set_title(v242group.title)
        #        print "group.start_pdf="+group.start_pdf
        #        print "group.npdf="+str(group.npdf)
        group.set_end_pdf(int(group.start_pdf)+int(group.npdf)-1)
        if "Centralscalevariation" in v242group.title :
            group.set_isscale(1)
            group.set_title("scale_variation")
            group.set_end_pdf(group.start_pdf)
            #            print "isscale"
        if "unknown" in v242group.title :
            group.set_title("PDF_"+str(group.start_pdf))
            
            
        groups.append(group)
        #print "##############"
        #print "title="+group.title
        #print "start_pdf"+str(group.start_pdf)
        #print "end_pdf"+str(group.end_pdf)
        
    return groups
#def make_template():
#    LO242=runLO242()
#    template=define_set(LO242.Group_list)
#    for group in template:
#        print "=============="
#        print "FROM : "+str(group.start_pdf)
#        if(not bool(group.isscale)): print "TO : " + str(int(group.start_pdf)+group.npdf-1)
#        else : print "SCALE VAR"
#    return template





    
    
    
def arrange_with_template(sample,template):
    weightlist=[]
    central_pdf=''
    for group_t in template:
        if bool(group_t.isscale): central_pdf=group_t.start_pdf
    for group_t in template:
        ##Let's do for scale variance case
#        if bool(group_t.isscale) :
            
        for group_s in sample.Group_list:
            for weight_s in group_s.weights:
                #print "##############"
                #print "weight_s.pdf="+weight_s.pdf
                #print "weight_s.id="+str(weight_s.id)
                #print "group_s.isscale="+str(group_s.isscale)
                #print "group_t.start_pdf="+group_t.start_pdf
                if ( "dyn_scale_choice" in weight_s.name ) :
                    continue
                if bool(group_s.isscale):###this is for NLO242
                    weight_s.pdf=central_pdf
                    #     if( (int(weight_s.pdf) == int(group_t.start_pdf) ) and bool(group_t.isscale) ) :
                    if ( bool(group_t.isscale) ):
                        info=weightinfo()
                        info.id=weight_s.id
                        info.pdf=weight_s.pdf
                        if( info.pdf==''): info.pdf=group_t.start_pdf
                        info.muF=weight_s.muF
                        info.muR=weight_s.muR
                        info.name="scale_variation_muR_"+str(weight_s.muR)+"_muF_"+str(weight_s.muF)
                        weightlist.append(info)
                        #                        print info.id
                elif( (int(weight_s.pdf) >= int(group_t.start_pdf)) and (int(weight_s.pdf) <= int(group_t.end_pdf)) ):
                    #                    print "pdf replica"+weight_s.name
                    info=weightinfo()
                    info.id=weight_s.id
                    info.pdf=weight_s.pdf
                    info.muF=weight_s.muF
                    if info.muF=='' : info.muF='1.0'
                    info.muR=weight_s.muR
                    if info.muR=='' : info.muR='1.0'
                    info.name=group_t.title
                    weightlist.append(info)
    print "#of weights="+str(len(weightlist))
    return weightlist


def check_list(list_):
    for weight in list_:
        print "["+weight.id+"] "+weight.name
                    
def print_list_cpp(list_,version):

    for weight in list_:
        print version+"Info.push_back(weightinfo("+'"'+str(weight.id)+'"'+", "+'"'+str(weight.pdf)+'"'+', "'+str(weight.name)+'", '+str(float(weight.muR))+", "+str(float(weight.muF))+"));"


if __name__ == "__main__":


    '''
    LO242=runLO242()
    template=define_set(LO242.Group_list)    
    list242=arrange_with_template(LO242,template)
    print "# of elements in v242="+str(len(list242))
    #    check_list(list242)
    print_list_cpp(list242,"LO242")
    LO26x=runLO26x()
    list26x=arrange_with_template(LO26x,template)
    print "# of elements in v26x="+str(len(list26x))
    '''
    #    check_list(list26x)
    #    print_list_cpp(list26x,"LO26x")

    NLO242=runNLO242()
    Nweight242=0
    for group in NLO242.Group_list:
#        print "##############"
 #       print "title="+group.title
  #      print "start_pdf="+group.weights[0].pdf
        #print "end_pdf="+group.end_pdf
   #     print "npdf="+str(len(group.weights))
        for weight in group.weights:
            #print "###################"
            #print "id="+str(weight.id)
            #print "pdf="+weight.pdf
            Nweight242+=1
    print "Nweight242="+str(Nweight242)

    NLO26x=runNLO26x()
    Nweight26x=0
    for group in NLO26x.Group_list:
        #        print "###############"
        #       print group.title
        for weight in group.weights:
            #print "id="+weight.id+"pdf="+weight.pdf
            Nweight26x+=1
    print "Nweight26x="+str(Nweight26x)
    #    template=define_set(NLO26x.Group_list)

    template=define_set(NLO242.Group_list,'306000')
    #for group in template:
        #print "##########"
        #print "title="+group.title
        #print "start_pdf="+str(group.start_pdf)
        #print "end_pdf="+str(group.end_pdf)
        #print "npdf="+str(group.npdf)
    list242=arrange_with_template(NLO242,template)
    #    for weight in list242:
    
    #       print "id="+weight.id+"name="+weight.name+" pdf="+weight.pdf
    list26x=arrange_with_template(NLO26x,template)
    #for weight in list26x:
    #    print weight.id
    #   list242=arrange_with_template(NLO242,template)
    #Group_list2=NLO26x.Group_list
    #for group in Group_list2:
    #    print group.title
#    print_list_cpp(list242,"NLO242")
 #   print_list_cpp(list26x,"NLO26x")

