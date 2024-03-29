#!/usr/bin/python3
import sys
import subprocess
import argparse
p = argparse.ArgumentParser(epilog="Created by Gareth Nisbet 2014",description="example: ggears -n 36 -p 20 -t internal -dia 44")
p.add_argument("-n", dest="num", help="Number of teeth (Required)")
p.add_argument('-p', dest='pressureangle', help = 'pressure Angle (Required)')
p.add_argument('-t', dest='type', help = 'internal or external (Required)')
p.add_argument('-dia', dest = 'diameter', help ='Diameter of hole for external gear or outer diameter of internal gear (Required)')
p.add_argument('-a', dest = 'addendum', help ='Manual override of addendum (normally 1/Pitch Diameter)')
p.add_argument('-d', dest = 'dedendum', help ='Manual override of dedendum (normally 1.25/Pitch Diameter)')
p.add_argument('-m', dest = 'module', help ='Module')
p.add_argument('-r', dest = 'resolution', help ='Sets the number of points in the involution')
p.add_argument('-f', dest = 'prefix', help ='Choose file prefix (example: -f name)')
p.add_argument('-tm', dest = 'tmat', help ='Transformation Matrix for Inkscape (default: 3.5434,0,0,3.5434,0,0 for mm scale)')
p.add_argument('-ft', dest = 'filetype', help ='(svg or dxf) Use this option if you want to bypass plotting or if Tkinter is not available')

args = p.parse_args()
try:
    import numpy as np
except:
    print('Numpy is not installed; try:\n sudo apt-get install numpy')
    exit()
try:
    try:
        from Tkinter import *
    except:
        from tkinter import *
except:
    print('Tkinter is not installed; try:\n sudo apt-get install Tkinter\nOr use -ft svg or -ft dxf to write files directly.')
    exit()
if not len(sys.argv) > 1:
    print('No arguments given. Type: ggears -h')
    exit()

numteeth = int(round(float(args.num)))
pressureangle=float(args.pressureangle)
base=float(args.diameter)
internal = args.type
if args.prefix == None:
    prefix = 'ggears-'
else:
    prefix = args.prefix
if args.resolution == None:
    resolution = 20
else:
    resolution=float(args.resolution)
outfile = args.prefix
if args.module == None:
    m=1
else:
    m=float(args.module)
    
pitchdia = numteeth*m
Pd = numteeth / pitchdia

if args.addendum == None:
    addendum = 1.0/Pd
else:
    addendum = float(args.addendum)
if args.dedendum == None:
    dedendum = 1.25/Pd
else:
    dedendum = float(args.dedendum)

od=pitchdia+(2.0*addendum)
ld=pitchdia-(2.0*addendum)
bd=pitchdia*np.cos(pressureangle*np.pi/180.0)

maxrange=np.sqrt(((od**2)/(bd**2))-1)
rootd=pitchdia-(2*dedendum)
if bd < rootd:
    minrange=np.sqrt(((rootd**2)/(bd**2))-1)
else:
    minrange=0

table = {'1 Root Diameter': rootd, '2 Base Diameter':  bd, '3 Pitch Diameter': pitchdia, \
         '4 Outside Diameter': od, '5 Addendum': addendum,'6 Dedendum': dedendum,'7 Pressure Angle': pressureangle,'8 Module': m,\
         '9 Resolution': resolution}
table2='1   No. Teeth: ' + str(numteeth) +'\n'+'2   Root Diameter: ' + str(rootd) +'\n'+ '3   Base Diameter: ' + str(bd) +'\n' + '4   Pitch Diameter: ' + str(pitchdia) +'\n' +'5   Outside Diameter:'+ str(od) +'\n' + '6   Addendum: '+ str(addendum)+'\n'+'7   Dedendum: '+str(dedendum) +'\n'+'8   Pressure Angle: '+str(pressureangle) +'\n'+'9   Module: '+ str(m) +'\n'+'10 Resolution: ' + str(resolution)
for param, val in sorted(table.items()):
    print('%-20s  %0.4f' % (param, val))

#===============================================================================
#                 Generalised rotation matrix
#===============================================================================

class rotxyz(object):
    """Example p = rotxyz(initial_vector, vectorrotateabout, angle)"""
    def __init__(self,u,angle):
        self.u = u
        self.angle = angle
        u=np.matrix(self.u)/np.linalg.norm(np.matrix(self.u))
        e11=u[0,0]**2+(1-u[0,0]**2)*np.cos(angle*np.pi/180.0)
        e12=u[0,0]*u[0,1]*(1-np.cos(angle*np.pi/180.0))-u[0,2]*np.sin(angle*np.pi/180.0)
        e13=u[0,0]*u[0,2]*(1-np.cos(angle*np.pi/180.0))+u[0,1]*np.sin(angle*np.pi/180.0)
        e21=u[0,0]*u[0,1]*(1-np.cos(angle*np.pi/180.0))+u[0,2]*np.sin(angle*np.pi/180.0)
        e22=u[0,1]**2+(1-u[0,1]**2)*np.cos(angle*np.pi/180.0)
        e23=u[0,1]*u[0,2]*(1-np.cos(angle*np.pi/180.0))-u[0,0]*np.sin(angle*np.pi/180.0)
        e31=u[0,0]*u[0,2]*(1-np.cos(angle*np.pi/180.0))-u[0,1]*np.sin(angle*np.pi/180.0)
        e32=u[0,1]*u[0,2]*(1-np.cos(angle*np.pi/180.0))+u[0,0]*np.sin(angle*np.pi/180.0)
        e33=u[0,2]**2+(1-u[0,2]**2)*np.cos(angle*np.pi/180.0)
        self.rotmat = np.matrix([[e11,e12,e13],[e21,e22,e23],[e31,e32,e33]])
    def rmat(self):
        return self.rotmat
    
#===============================================================================
#                         involute creator
#===============================================================================
def involute(dedendum,bd,pd,datarange):
    pmax=np.sqrt(((pd**2)/(bd**2))-1)
    xp1=(bd/2.0)*(np.cos(pmax)+(pmax)*np.sin(pmax))
    yp1=(bd/2.0)*(np.sin(pmax)-(pmax)*np.cos(pmax))
    xp2=-(bd/2.0)*(np.cos(pmax-(np.pi/2))+(pmax)*np.sin(pmax-(np.pi/2)))
    yp2=-(bd/2.0)*(np.sin(pmax-(np.pi/2))-(pmax)*np.cos(pmax-(np.pi/2)))
    xs=np.array([(bd/2.0)*(np.cos(datarange)+(datarange)*np.sin(datarange))])-xp1
    ys=np.array([(bd/2.0)*(np.sin(datarange)-(datarange)*np.cos(datarange))])-yp1
    zs=np.array([np.zeros(len(datarange))])
    ys2=np.array([-(bd/2.0)*(np.cos(datarange-(np.pi/2))+(datarange)*np.sin(datarange-(np.pi/2)))])-xp2
    xs2=np.array([-(bd/2.0)*(np.sin(datarange-(np.pi/2))-(datarange)*np.cos(datarange-(np.pi/2)))])-yp2
    zs2=np.array([-np.zeros(len(datarange))])
    bv=np.array([[((dedendum/2.0)-xp1)]])
    xs=np.hstack((bv,xs))
    ys=np.array([np.hstack((ys[0][0],ys[0]))])
    zs=np.array([np.hstack((zs[0][0],zs[0]))])
    xs2=np.hstack((bv,xs2))
    ys2=np.array([np.hstack((ys2[0][0],ys2[0]))])
    zs2=np.array([np.hstack((zs2[0][0],zs2[0]))])
    return np.concatenate((xs.T,ys.T,zs.T,xs2.T,ys2.T,zs2.T),1)

#===============================================================================
#             Gear Creator
#===============================================================================
def gears(numteeth, rootd, pitchdia,w1,invol,gtype):
    invols=np.array([[np.NAN,np.NAN,np.NAN]])
    pi=np.pi
    if gtype == 'internal' or gtype == '1':
        w1=-w1
        for i1 in np.arange(0,-360.0,-360./numteeth):
            bwallr=(rootd/2.0)
            pwallr=(pitchdia/2.0)
            r1=rotxyz([0,0,1],i1+w1+180).rmat()
            r2=rotxyz([0,0,1],i1+180).rmat()
            invol1=np.flipud(invol[:,[0,1,2]]*r1)
            invol2=invol[:,[3,4,5]]*r2
            basewall1=np.array([[(bwallr)*np.cos(-i1*pi/180),(bwallr)*np.sin(-i1*pi/180),0]])
            basewall2=np.array([[(bwallr)*np.cos(-(i1+w1)*pi/180),(bwallr)*np.sin(-(i1+w1)*pi/180),0]])
            pitchwall1=np.array([[(pwallr)*np.cos(-i1*pi/180),(pwallr)*np.sin(-i1*pi/180),0]])
            pitchwall2=np.array([[(pwallr)*np.cos(-(i1+w1)*pi/180),(pwallr)*np.sin(-(i1+w1)*pi/180),0]])
            invol1=np.array(invol1+pitchwall1)
            invol2=np.array(invol2+pitchwall2)
            basewall1=np.array([[(rootd/2.0)*np.cos(-i1*pi/180),(rootd/2.0)*np.sin(-i1*pi/180),0]])
            basewall2=np.array([[(rootd/2.0)*np.cos(-(i1+w1)*pi/180),(rootd/2.0)*np.sin(-(i1+w1)*pi/180),0]])
            invol1=np.flipud(invol1)
            invol2=np.flipud(invol2)
            invols1=np.concatenate((invol1,invol2),0)
            invols=np.concatenate((invols1,invols),0)
    else:
        for i1 in np.arange(0,-360.0,-360./numteeth):
            bwallr=(rootd/2.0)
            pwallr=(pitchdia/2.0)
            r1=rotxyz([0,0,1],i1).rmat()
            r2=rotxyz([0,0,1],i1+w1).rmat()
            invol1=np.flipud(invol[:,[0,1,2]]*r1)
            invol2=invol[:,[3,4,5]]*r2
            basewall1=np.array([[(bwallr)*np.cos(-i1*pi/180),(bwallr)*np.sin(-i1*pi/180),0]])
            basewall2=np.array([[(bwallr)*np.cos(-(i1+w1)*pi/180),(bwallr)*np.sin(-(i1+w1)*pi/180),0]])
            pitchwall1=np.array([[(pwallr)*np.cos(-i1*pi/180),(pwallr)*np.sin(-i1*pi/180),0]])
            pitchwall2=np.array([[(pwallr)*np.cos(-(i1+w1)*pi/180),(pwallr)*np.sin(-(i1+w1)*pi/180),0]])
            invol1=np.array(invol1+pitchwall1)
            invol2=np.array(invol2+pitchwall2)
            basewall1=np.array([[(rootd/2.0)*np.cos(-i1*pi/180),(rootd/2.0)*np.sin(-i1*pi/180),0]])
            basewall2=np.array([[(rootd/2.0)*np.cos(-(i1+w1)*pi/180),(rootd/2.0)*np.sin(-(i1+w1)*pi/180),0]])
            invols1=np.concatenate((invol2,invol1),0)
            invols=np.concatenate((invols1,invols),0)
    invols=invols[~np.isnan(invols[:,0]),:]
    invols=np.concatenate((invols,invols[[0,1],:]),0)
    return invols
#===============================================================================
#                 DXF Creator
#===============================================================================
def createdxf():
    try:
        outfiledxf=outfile2+'.dxf'
        print('DXF created ' + outfile2 +'.dxf')
    except:
        outfile = 'internal-'+str(numteeth)+'teeth'+'M'+str(m)
        outfiledxf=outfile+'.dxf'
        print('DXF created ' + outfile +'.dxf')
    f = open(outfiledxf, "w")
    f.write("  0\nSECTION\n")
    f.write("  2\nENTITIES\n")
    f.write("  0\nLWPOLYLINE\n")
    f.write("  100\nAcDbPolyline\n")
    f.write("  90\n")
    f.write("  "+str(len(invols)-1)+"\n")
    f.write("  70\n0\n")
    for i1 in range(0,len(invols)):
        x1 = invols[i1,0]
        y1 = invols[i1,1]
        f.write(" 10\n{0}\n".format(x1))
        f.write(" 20\n{0}\n".format(y1))
        f.write(" 30\n0\n")
    f.write(" 10\n{0}\n".format(invols[0,0]))
    f.write(" 20\n{0}\n".format(invols[0,1]))
    f.write(" 30\n0\n")

    circvals=np.linspace(-np.pi,np.pi,120)
    circ=np.array([(base/2.0)*np.cos(circvals),(base/2.0)*np.sin(circvals)])
    f.write("  0\nLWPOLYLINE\n")
    f.write("  100\nAcDbPolyline\n")
    f.write("  90\n")
    f.write("  "+str(circ.shape[1])+"\n")
    f.write("  70\n0\n")
    for i1 in list(range(0,circ.shape[1])):
        cx1 = circ[0,i1]
        cy1 = circ[1,i1]
        f.write(" 10\n{0}\n".format(cx1))
        f.write(" 20\n{0}\n".format(cy1))
        f.write(" 30\n0\n")
    f.write("  0\nENDSEC\n")    
    f.write("  0\nEOF\n")
    f.close()
#===============================================================================
#             SVG Creator
#===============================================================================
def createsvg():
    try:
        outfilesvg=outfile2+'.svg'
        print('SVG created ' + outfile2 +'.svg')
    except:
        outfile = 'internal-'+str(numteeth)+'teeth'+'M'+str(m)
        outfilesvg=outfile+'.svg'
        print('SVG created ' + outfile +'.svg')
    f = open(outfilesvg, "w")
    f.write('<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n')
    f.write('<!-- Created with ggears by Gareth Nisbet) -->\n')
    f.write('<svg')
    f.write(' docname="ggears.svg">\n')
    f.write('   <g\n')
    if args.tmat == None:
        f.write('     transform="matrix(3.5434,0,0,3.5434,0,0)"\n')
    else:
        f.write('     transform="matrix('+str(args.tmat) +')"\n')
    f.write('     style="font-size:1px;stroke-width:0.05">\n')
    if internal=='internal':
    	f.write('<circle r="' +str(base/2.0) +'" stroke="black" stroke-width="0.05" fill="none" />\n')
    f.write('    <path\n')
    f.write('       id="gear"\n')
    f.write('       d="m ')
    for i in list(range(1,len(invols)-1)):
        x1 = invols[i-1,0]
        y1 = invols[i-1,1]
        f.write(str(x1)+','+str(y1)+' L ')
    f.write(str(invols[0,0])+','+str(invols[0,1]))
    f.write('z"\n')
    f.write('       style="fill:none;stroke:#000000" />\n')
    if internal=='external':
    	f.write('<circle r="' +str(base/2.0) +'" stroke="black" stroke-width="0.05" fill="none" />\n')
    f.write('  </g>\n')
    f.write('</svg>')
    f.close()
    
#===============================================================================
#         GUI Commands
#===============================================================================
def inkopen():
    createsvg()
    try:
        subprocess.Popen(['inkscape', outfile2 + '.svg'], stdout=subprocess.PIPE)
    except:
        outfile = 'internal-'+str(numteeth)+'teeth'+'M'+str(m)
        subprocess.Popen(['inkscape', outfile + '.svg'], stdout=subprocess.PIPE)
def freeCADopen():
    createdxf()
    try:
        subprocess.Popen(['freecad', outfile2 +'.dxf'], stdout=subprocess.PIPE)
    except:
        outfile = 'internal-'+str(numteeth)+'teeth'+'M'+str(m)
        subprocess.Popen(['freecad', outfile +'.dxf'], stdout=subprocess.PIPE)
def exitfunc():
    exit()
#===============================================================================
#         Main Gui
#===============================================================================
class Gears( Frame ):
    def __init__( self ):
        Frame.__init__( self )
        self.grid( columnspan=20, rowspan=20 )
        self.master.title( 'GGears v 0.1' )
        self.master.geometry( "770x600" )
        self.control = Scale( self, from_ = -3, to = 100, resolution=1 ,orient = HORIZONTAL, command = self.scaleGear )
        self.control.grid(row=0, column=0, columnspan=10, rowspan=20,sticky=W+E+N+S )
        self.control.set( 0 )
        self.canvasx = 550
        self.canvasy = 550
        self.display = Canvas( self, bg = "black", width=self.canvasx, height=self.canvasy)
        self.display.grid(columnspan=10, rowspan=20 )
#===============================================================================
#         Buttons
#===============================================================================
        self.button = Button(text="Incscape", width=10, fg="black", command = inkopen)
        self.button.grid(row=3, column=22)
        self.button = Button(text="freeCAD", width=10,fg="black", command = freeCADopen)
        self.button.grid(row=4, column=22)
        self.button = Button(text="Create DXF", width=10,fg="black", command = createdxf)
        self.button.grid(row=3, column=23)
        self.button = Button(text="Create SVG",width=10, fg="black", command = createsvg)
        self.button.grid(row=4, column=23)
#         self.button = Button(text="Exit", width=10,fg="red",command = exitfunc)
#         self.button.grid(row=18, column=22)
        
        self.label = Label(self.master, text= "No. Teeth")
        self.label.grid(row=7, column=22)
        self.entrytext = StringVar()
        Entry(self.master, textvariable=self.entrytext,width=13).grid(row=7, column=23)
        self.buttontext = StringVar()
        self.label = Label(self.master, text= "Pressure Angle")
        self.label.grid(row=8, column=22)
        self.entrytext2 = StringVar()
        Entry(self.master, textvariable=self.entrytext2,width=13).grid(row=8, column=23)
        self.buttontext = StringVar()
        
        self.label = Label(self.master, text= "Adendum")
        self.label.grid(row=9, column=22)
        self.entrytext3 = StringVar()
        Entry(self.master, textvariable=self.entrytext3,width=13).grid(row=9, column=23)
        self.buttontext = StringVar()
        
        self.label = Label(self.master, text= "Dedendum")
        self.label.grid(row=10, column=22)
        self.entrytext4 = StringVar()
        Entry(self.master, textvariable=self.entrytext4,width=13).grid(row=10, column=23)
        self.buttontext = StringVar()
        
        self.label = Label(self.master, text= "Hole Size")
        self.label.grid(row=11, column=22)
        self.entrytext5 = StringVar()
        Entry(self.master, textvariable=self.entrytext5,width=13).grid(row=11, column=23)
        self.buttontext = StringVar()
        
        self.label = Label(self.master, text= "Resolution")
        self.label.grid(row=12, column=22)
        self.entrytext6 = StringVar()
        Entry(self.master, textvariable=self.entrytext6,width=13).grid(row=12, column=23)
        self.buttontext = StringVar()
        
        self.label = Label(self.master, text= "Module")
        self.label.grid(row=13, column=22)
        self.entrytext7 = StringVar()
        Entry(self.master, textvariable=self.entrytext7,width=13).grid(row=13, column=23)
        self.buttontext = StringVar()
        
        self.var1 = IntVar()
        self.var1=StringVar()
        Checkbutton(self.master, width=10, text="Internal", variable=self.var1).grid(row=14, column=23)
        
        self.buttontext.set("Update")
        Button(self.master, textvariable=self.buttontext,width=10, fg="red", command=self.onClick).grid(row=15, column=23)

    def onClick(self):
        if self.entrytext.get() == '':
            numteeth2 = numteeth
        else:
            numteeth2 = float(self.entrytext.get())
        if self.entrytext2.get() == '':
            pressureangle2=pressureangle
        else:
            pressureangle2 = float(self.entrytext2.get())
        global base
        if self.entrytext5.get() == '':
            base2 = base
        else:
            base2 = float(self.entrytext5.get())
        global resolution2
        if self.entrytext6.get() == '':
            resolution2 = resolution
        else:
            resolution2 = float(self.entrytext6.get())
        
        global m2
        if self.entrytext7.get() == '':
            m2 = m
        else:
            m2 = float(self.entrytext7.get())

        if self.entrytext3.get() == '':
            pitchdia2 = numteeth2*m2
            Pd2 = numteeth2 / pitchdia2
            addendum2 = 1.0/Pd2
        else:
            addendum2 = float(self.entrytext3.get())

        if self.entrytext4.get() == '':
            pitchdia2 = numteeth2*m2
            Pd2 = numteeth2 / pitchdia2
            dedendum2 = 1.25/Pd2

        else:
            dedendum2 = float(self.entrytext4.get())
        
        pitchdia2 = numteeth2*m2
        Pd2 = numteeth2 / pitchdia2
        bd2=pitchdia2*np.cos(pressureangle2*np.pi/180.0)
        od2=pitchdia2+(2.0*addendum2)
        w2=-180.0/numteeth2
        rootd2=pitchdia2-(2*dedendum2)
        if bd2 < rootd2:
            minrange=np.sqrt(((rootd2**2)/(bd2**2))-1)
        else:
            minrange=0
        maxrange=np.sqrt(((od2**2)/(bd2**2))-1)
        global gtype
        gtype=self.var1.get()
        self.display.delete( "gear" )
        canvasx=self.canvasx
        canvasy=self.canvasy
        
        base=base2       
        global invols

        datarange2=np.linspace(minrange,maxrange,int(resolution2))
        invol2=involute(rootd2,bd2,pitchdia2,datarange2)
        invols=gears(numteeth2, rootd2, pitchdia2,w2,invol2,gtype)
        for i in range(1,len(invols)-1):
            x1 = invols[i-1,0]
            y1 = invols[i-1,1]
            x2 = invols[i,0]
            y2 = invols[i,1]
            self.display.create_line( canvasx/2.0 + x1*scale, canvasy/2.0 - y1*scale, canvasx/2.0 + x2*scale, canvasy/2.0 - y2*scale, fill="#acff77", tag = 'gear')
            x = base*scale/2.0
        self.display.create_oval( canvasx/2.0 - x, canvasy/2.0 - x, canvasx/2.0 + x, canvasy/2.0 + x, outline="#acff77", tag='gear')
        canvas_id=self.display.create_text(10, 10, anchor="nw",fill='white')
        global table2
        table2='1   No. Teeth: ' + str(numteeth2) + '\n'+ '2   Root Diameter: ' + str(rootd2) +'\n'+ '3   Base Diameter: ' + str(bd2) +'\n' + '4   Pitch Diameter: ' + str(pitchdia2) +'\n' +'5   Outside Diameter:'+ str(od2) +'\n' + '6   Addendum: '+ str(addendum2)+'\n'+'7   Dedendum: '+str(dedendum2) +'\n'+'8   Pressure Angle: '+str(pressureangle2) +'\n'+'9   Module: '+ str(m2) +'\n'+'10 Resolution: ' + str(resolution2)
        self.display.itemconfig(canvas_id, text=table2, tag='gear')
        global outfile2
        if gtype == 'internal' or gtype == '1':
            outfile2 = 'internal-'+str(numteeth2)+'teeth'+'M'+str(m)
        else:
            outfile2 = prefix + 'External-'+str(numteeth2)+'teeth'+'M'+str(m)
#         self.master.destroy()
#         print(input)
#     def button_click(self, e):
#         pass

    def scaleGear( self, scaleValue):
        global scale
        scale = int( scaleValue ) + 3.07
        self.display.delete( "gear" )
        canvasx=self.canvasx
        canvasy=self.canvasy
        for i in list(range(1,len(invols)-1)):
            x1 = invols[i-1,0]
            y1 = invols[i-1,1]
            x2 = invols[i,0]
            y2 = invols[i,1]
            self.display.create_line( canvasx/2.0 + x1*scale, canvasy/2.0 - y1*scale, canvasx/2.0 + x2*scale, canvasy/2.0 - y2*scale, fill="#acff77", tag = 'gear',width=1)
            x = base*scale/2.0
        self.display.create_oval( canvasx/2.0 - x, canvasy/2.0 - x, canvasx/2.0 + x, canvasy/2.0 + x, outline="#acff77", tag='gear',width=1)
        canvas_id=self.display.create_text(10, 10, anchor="nw",fill='white')
        self.display.itemconfig(canvas_id, text=table2, tag='gear')

#===============================================================================
#             Run the Program
#===============================================================================
datarange=np.linspace(minrange,maxrange,resolution)
invol=involute(rootd,bd,pitchdia,datarange)
w1=-180.0/numteeth
invols=gears(numteeth, rootd, pitchdia,w1,invol,args.type)

if args.filetype != None:
    if args.type == 'internal':
        outfile = prefix + '-internal-'+str(numteeth)+'teeth'+'M'+str(m)
    else:
        outfile = prefix + '-external-'+str(numteeth)+'teeth'+'M'+str(m)
    if args.filetype == 'svg':
        outfile = outfile+'.svg'
        createsvg()
        exit()
    elif args.filetype == 'dxf':
        outfile = outfile+'.dxf'
        createdxf()
        exit()
def main():
   Gears().mainloop()   

if __name__ == "__main__":
   main()
