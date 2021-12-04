##########################
#   SISO TOOL in Python  #
# JP Nouaille - dec 2019 #
##########################
# needs: matplotlib, control, scipy, pylab, numpy, math
from matplotlib.pyplot import *
import control 
from control.matlab import *
from math import *
from scipy.interpolate import interp1d
from pylab import *
from matplotlib.widgets import Slider
import numpy as np
import matplotlib
#matplotlib.interactive(True)

#backend='Qt5Agg'
#backend='Qt4Agg'
#backend='TkAgg'
#backend='WXAgg'
#matplotlib.rcParams['backend'] = backend

#switch_backend(backend)

def step_info(t,yout):
    """ Return info on a step response yout as a function    
        of time t:    
        Overshoot OS (%)    
        Rise time at 66 % of the steady state response    
        Settling time at 5 % of the steady state response """

   
    OS=(yout.max()/yout[-1]-1.0)*100.0
    Tr=t[next((i for i in range(0,len(yout)-1) if yout[i]>yout[-1]*.66),len(yout))-1]-t[0]
    Ts=t[next((len(yout)-i for i in range(2,len(yout)-1) if abs(yout[-i]/yout[-1])>1.05 or abs(yout[-i]/yout[-1])<0.95 ),len(yout))-1]-t[0]
    # Console output
    #print( "OS: %f%s"%((yout.max()/yout[-1]-1)*100,'%'))
    #print( "Tr: %fs"%(t[next(i for i in range(0,len(yout)-1) if yout[i]>yout[-1]*.66)]-t[0]))
    #print( "Ts: %fs"%(t[next(len(yout)-i for i in range(2,len(yout)-1) if abs(yout[-i]/yout[-1])>1.05 or abs(yout[-i]/yout[-1])<0.95 )]-t[0]))
    return (OS,Tr,Ts)


def damp(sys):
    """ Return, for each pole of the system  
        The pole ri
        The real part a    
        The imaginary part b    
        The damping ratio    
        The pulsation (rad/s) 
        A string for output with a,b,xi,w """
    roots=control.pole(sys)
    ri=[]
    a=[]
    b=[]
    w=[]
    xi=[]
    st=[]
    for i in range(0,roots.size):
      #print('i=',i)
      ri.append(roots[i])
      a.append(roots[i].real)
      b.append(roots[i].imag)
      w.append(math.sqrt(a[i]**2+b[i]**2))
      xi.append(-a[i]/w[i])
      if b[i]>0:
        signb='+'
      else:
        signb='-' 
      #print(repr(a[i])+signb+'j'+repr(math.fabs(b[i]))+'  xi='+repr(xi[i])+' w='+repr(w[i])+' rad/s')
      st.append('%.3f'%(a[i])+signb+'j'+'%.3f'%(math.fabs(b[i]))+'  xi='+'%.3f'%(xi[i])+'  w='+'%.3f'%(w[i])+' rad/s')
      #print(st[i])
    return (ri,a,b,xi,w,st)



class dragGUI:
    def __init__(self,sys=-control.tf([-13.76,-8.466],[1.,1.569,13.86]),mingain=0.0001,maxgain=40.0,kdefault=1.0):
        self.fig1=figure(1)
        self.fig1.clf()
        self.ax1=self.fig1.add_subplot(222)
        subplots_adjust(left=0.1, bottom=0.20, right=0.95,top=0.95)

        #sysqdelta=control.tf([-13.76,-8.466],[1.,1.569,13.86])
        #self.sys=-sysqdelta
        #sysgammabo=control.tf([ -0.05088635,  -0.04604676,  11.02724702],[  1.00000000e+00,   1.13128903e+01,   7.67080677e+01, -4.03005685e-15])
        #self.sys=sysgammabo
        # pour z
        #sysz=control.tf(array([  -353.8228632 ,   -320.17225056,  76674.63655141]),array([  1.00000000e+00,   1.04406983e+01,   7.59188263e+01,
        #   1.89007014e+02,  -4.36558782e-11]))
        #self.sys=sysz
        self.sys=sys


        #mingain=0.0001
        #maxgain=40
        if kdefault==1:
            self.gain=(mingain+maxgain)/2.0
        else:
            self.gain=kdefault

        (roots,gains) =control.rlocus(self.sys,logspace(log10(mingain),log10(maxgain),3000),Plot=False)
        self.roots=roots
        
        self.rootlocus=self.ax1.plot(roots.real,roots.imag)

        polesloc=self.ax1.plot(roots.real,roots.imag)
        rootsol=pole(self.sys)
        zerool=zero(self.sys)
        self.ax1.plot(rootsol.real,rootsol.imag,'bx',markersize=10) 
        self.ax1.plot(zerool.real,zerool.imag,'bo',markersize=10) 


        clsys=control.feedback(series(tf([self.gain],[1]),self.sys),tf([1],[1]))
        rootsk0=control.pole(clsys)
        self.polesloc, =self.ax1.plot(rootsk0.real,rootsk0.imag,'r*',markersize=10)

        aa=max(sqrt(rootsk0.real**2+rootsk0.imag**2))
        # xi=0.7 axes
        self.ax1.plot([0,-aa],[0,aa],'k--')
        self.ax1.plot([0,-aa],[0,-aa],'k--')

        #self.ax1.title='Root locus'
        self.ax1.title.set_text('Root locus')  
        self.ax1.grid(True)
        #self.polesloc,=self.ax1.plot([1],[1],'r*',lw=2,markersize=10)

        (y,t)=control.matlab.step(clsys)
        # Time response on figure 1, subplot 2
        self.ax2=self.fig1.add_subplot(224)
        self.l2, =self.ax2.plot(t,y,'b',lw=2)
        self.ax2.set_xlabel('Time (s)')
        self.ax2.set_ylabel('y(t)')
        self.ax2.set_autoscale_on
        self.ax2.grid(True)

        # check system sys stability
        stable=True
        for i in range(0,rootsk0.size):
          if (rootsk0[i].real)>=0:
            stable=False
        if stable==True:
          (OS,Tr,Ts)=step_info(t,y)
          self.l5=figtext(0.1,0.1,'OS='+'{:.3f}'.format(OS)+' %'+'   tr5%='+'{:.3f}'.format(Ts)+' s'+'   Gain='+'{:.5f}'.format(self.gain))
          yy=interp1d(t,y)
          self.l21, =self.ax2.plot(Ts,yy(Ts),'r*')
        else:
          self.l5=figtext(0.1,0.1,'OS=N/A %   tr2%=N/A s'\
                       +'   Gain='+'{:.5f}'.format(self.gain))
          self.l21, =self.ax2.plot(t[-1],y[-1],'r*')
        # time response
        self.l23,=self.ax2.plot([t[0],t[-1]],[y[-1]*0.95,y[-1]*0.95],'k--')
        self.l24,=self.ax2.plot([t[0],t[-1]],[y[-1]*1.05,y[-1]*1.05],'k--')

        # open loop bode
        self.ax3=self.fig1.add_subplot(321)
        #self.fig2,(self.ax3,self.ax4) = subplots(2)
        [Gaindb,phasedeg,omega]=bode(self.gain*self.sys,Plot=False)
        self.l3,=self.ax3.loglog(omega,Gaindb,lw=2)
        #ax3.set_xlabel('Freq (rad/s)')
        self.ax3.set_ylabel('Gain (dB)')
        self.ax3.grid(True,which='both')
        self.ax3.set_title('Open loop Bode')
        
        
        self.ax4=self.fig1.add_subplot(323)
        self.ax4.set_xlabel('Freq (rad/s)')
        self.ax4.set_ylabel('Phase (deg)')
        self.ax4.grid(True,which='both')
        
        self.l4,=self.ax4.semilogx(omega,phasedeg,lw=2)

        # margin
        gm, pm, wg, wp = margin(self.gain*self.sys)
        self.l6=figtext(0.1,0.05,'GM='+'{:.2f}'.format(gm)+' dB   PM='+'{:.3f}'.format(pm)+' deg') 

        # closed loop poles characteristics
        (ri,a,b,xi,w,st)=damp(clsys)

        self.l7=[]
        for i in range(0,len(st)):
            self.l7p=figtext(0.02,0.3-i*0.035,st[i])
            self.l7.append(self.l7p)


        self.pl=draggableGain(self,self.ax1,roots,gains,self.gain)
        self.pl.connect()

    def update(self,x,y,val):
        """ Scale callback : update the 3 windows of the HMI """
        z=x+1j*y
        #print('z=%f'% real(z),'+(%f)j'% imag(z))
        roots2=self.roots-z
        #print('size(roots2)=(%i,%i)'%roots2.shape)
        noroots=roots2.shape[1]
        closestRoots=argmin(abs(roots2))

        i=int(closestRoots/noroots)
        j=int(closestRoots%noroots)
        #print('i=%d'%i,'j=%d'%j)
        #print('roots2[%d'%i,'%d'%j,']={:f}'.format(roots2[i,j]))
        deltax=abs(roots2[i,j])

        # get gain
        self.k = val
        clsys=control.feedback(series(tf([self.k],1),self.sys),tf([1],[1]))

        # closed loop poles
        rootsk=control.pole(clsys)
        #self.polesloc.set_xdata(x)
        #self.polesloc.set_ydata(y)
        self.polesloc.set_xdata(rootsk.real)
        self.polesloc.set_ydata(rootsk.imag)

        (y,t)=control.matlab.step(clsys)
        
        # check system sys stability
        stable=True
        for i in range(0,rootsk.size):
          if (rootsk[i].real)>=0:
            stable=False
        if stable==True:
          # step response characteristics (overshoot, settling time)
          (OS,Tr,Ts)=step_info(t,y)
          self.l5.set_text('OS='+'{:.3f}'.format(OS)+' %'+'   tr5%='+'{:.3f}'.format(Ts)+' s'\
                          +'   Gain='+'{:.5f}'.format(self.k))
          yy=interp1d(t,y)
          self.l21.set_xdata(Ts)
          self.l21.set_ydata(yy(Ts)) 
          self.l23.set_xdata([t[0],t[-1]])
          self.l23.set_ydata([0.95*y[-1],0.95*y[-1]])
          self.l24.set_xdata([t[0],t[-1]])
          self.l24.set_ydata([1.05*y[-1],1.05*y[-1]])
        else:
          self.l5.set_text('OS= N/A   tr5%= N/A s'+'   Gain='+'{:.5f}'.format(self.k))

        self.l2.set_xdata(t)
        self.l2.set_ydata(y)
        self.ax2.set_xlim([min(t),max(t)])
        self.ax2.set_ylim([min(y),max(y)*1.1])

        # poles characterisitics
        (ri,a,b,xi,w,st)=damp(clsys)

        # margin
        gm, pm, wg, wp = margin(self.k*self.sys)
        # console output
        #print('gm=',gm,' dB   pm=',pm,' deg')
    
        # open loop bode 
        [GaindB,phasedeg,omega]=bode(self.k*self.sys,Plot=False)


        # for change of the gain, the phase is unchanged
        self.l4.set_xdata(omega)
        self.l4.set_ydata(phasedeg)
        self.ax4.set_xlim([min(omega),max(omega)])
        self.ax4.set_ylim([min(phasedeg),max(phasedeg)])

        self.l6.set_text('GM='+'{:.3f}'.format(gm)+' dB   PM='+'{:.3f}'.format(pm)+' deg')
        for i in range(0,len(st)):
            self.l7[i]._text=st[i]

        # open loop bode 
        [GaindB,phasedeg,omega]=bode(self.k*self.sys,Plot=False)
        #l3.set_xdata(omega)
        #GaindB=Gaindb+20*log10(abs(k))-20*log10(abs(self.k0))
        self.l3.set_ydata(GaindB)
        self.ax3.set_xlim([min(omega),max(omega)])
        self.ax3.set_ylim([min(GaindB),max(GaindB)])
        
        # for change of the gain, the phase is unchanged
        self.l4.set_xdata(omega)
        self.l4.set_ydata(phasedeg)
        self.ax4.set_xlim([min(omega),max(omega)])
        self.ax4.set_ylim([min(phasedeg),max(phasedeg)])

    
        self.fig1.canvas.draw_idle()
 

class draggableGain:
    def __init__(self, GUI, polesloc, roots,gains, k):
        """ Init draggagle gain """
        self.GUI=GUI
        self.polesloc = polesloc
        self.press = None
        self.dx=0
        self.dy=0
        self.roots = roots
        self.gains = gains
        self.dims = roots.shape
        self.k=k


    def connect(self):
        """ connect to all the events we need """
        self.cidpress = self.polesloc.figure.canvas.mpl_connect(
            'button_press_event', self.on_press)
        self.cidrelease = self.polesloc.figure.canvas.mpl_connect(
            'button_release_event', self.on_release)
        self.cidmotion = self.polesloc.figure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)


    def on_press(self, event):
        """ Action when a mouse button is pressed. """
        toolbar = plt.get_current_fig_manager().toolbar
        if toolbar.mode!='': return
        'on button press we will see if the mouse is over us and store some data'
        #print('mouse left click')
        if event.inaxes != self.polesloc.axes: return
        if self.press == None :
            self.x0=event.xdata
            self.y0=event.ydata
            self.dx=0
            self.dy=0
        self.press = self.x0, self.y0, event.xdata, event.ydata
        z=event.xdata+1j*event.ydata
        #print('z=%f'% real(z),'+(%f)j'% imag(z))
        roots2=self.roots-z
        #print('size(roots2)=(%i,%i)'%roots2.shape)
        noroots=roots2.shape[1]
        closestRoots=argmin(abs(roots2))
        i=int(closestRoots/noroots)
        j=int(closestRoots%noroots)
        #print('i=%d'%i,'j=%d'%j)
        #print('roots2[%d'%i,'%d'%j,']={:f}'.format(roots2[i,j]))
        deltax=abs(roots2[i,j])
        if i<1:
          self.k=self.gains[0]
        elif i>=roots2.shape[0]:
          self.k=self.gains[-1]
        else:
          self.k=self.gains[i];
        self.GUI.update(self.x0+self.dx,self.y0+self.dy,self.k)

        

    def on_motion(self, event):
        'on motion we will move the polesloc if the mouse is over us'
        if self.press is None: return
        if event.inaxes != self.polesloc.axes: return
        self.x0, self.y0, xpress, ypress = self.press
        self.dx = event.xdata - self.x0
        self.dy = event.ydata - self.y0
        z=event.xdata+1j*event.ydata
        #print('z=%f'% real(z),'+(%f)j'% imag(z))
        roots2=self.roots-z
        #print('size(roots2)=(%i,%i)'%roots2.shape)
        noroots=roots2.shape[1]
        closestRoots=argmin(abs(roots2))
        i=int(closestRoots/noroots)
        j=int(closestRoots%noroots)
        #print('i=%d'%i,'j=%d'%j)
        #print('roots2[%d'%i,'%d'%j,']={:f}'.format(roots2[i,j]))
        deltax=abs(roots2[i,j])
        if i<1:
          self.k=self.gains[0]
        elif i>=roots2.shape[0]:
          self.k=self.gains[-1]
        else:
          self.k=self.gains[i];
        self.GUI.update(self.x0+self.dx,self.y0+self.dy,self.k)
        #print('gain = %f'%self.k)


    def on_release(self, event):
        'on release we reset the press data'
        if event.inaxes != self.polesloc.axes: return

        self.press = None
        self.x0=event.xdata
        self.y0=event.ydata
        self.dx=0
        self.dy=0

        self.GUI.update(self.x0+self.dx,self.y0+self.dy,self.k)


    def disconnect(self):
        'disconnect all the stored connection ids'
        self.polesloc.figure.canvas.mpl_disconnect(self.cidpress)
        self.polesloc.figure.canvas.mpl_disconnect(self.cidrelease)
        self.polesloc.figure.canvas.mpl_disconnect(self.cidmotion)
   
def sisotool(sys,kmin=0.0001,kmax=40.0,kdefault=1.0):
   sisotool.GUI=dragGUI(sys,kmin,kmax,kdefault)
   show()

if __name__ == '__main__': 
     GUI=dragGUI()
     show() 
