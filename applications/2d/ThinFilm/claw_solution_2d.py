import numpy
import os


class ClawSolution:

    def __init__(self, directory='_output', frame_number=0):

        self.frame_number = frame_number
        
        if directory[0] == '/':
            self.directory = directory
        else: 
            self.directory = os.getcwd() + '/' + directory


        #---- Read parameters from time file ----
        
        # time_file_name = self.timeFileName(frame_number)
        time_file = self.timeFile(frame_number)
        
        self.t    = float( time_file.readline().split()[0] )
        self.meqn = int( time_file.readline().split()[0] )
        time_file.readline()
        time_file.readline()
        self.dim  = int( time_file.readline().split()[0] )

        time_file.close()
        
        
        #---- Read parameters from data file ----

        # data_file_name = self.dataFileName(frame_number)
        data_file = self.dataFile(frame_number)
        
        data_file.readline()
        data_file.readline()
        self.mx    = int( data_file.readline().split()[0] )
        self.my    = int( data_file.readline().split()[0] )
        self.x_low = float( data_file.readline().split()[0] )
        self.y_low = float( data_file.readline().split()[0] )
        self.dx    = float( data_file.readline().split()[0] )
        self.dy    = float( data_file.readline().split()[0] )
        data_file.readline()
        
        self.q    = numpy.zeros( (self.mx, self.my, self.meqn) )

        for iy in range(self.my):
            
            for ix in range(self.mx):
                self.q[ix,iy,:] = numpy.array( map( float, data_file.readline().split() ) )

            data_file.readline()
                            
        data_file.close()
        
        
        
        #---- Other parameters ----
        
        self.x_high = self.x_low + self.mx * self.dx
        self.y_high = self.y_low + self.my * self.dy


    def frameString(self, frame_number):
        return '%04i' % frame_number


    def dataFile(self, frame_number):
        data_file_name = self.directory + '/fort.q' + self.frameString(frame_number)
        return open(data_file_name, 'r')
        

    def timeFile(self, frame_number):
        time_file_name = self.directory + '/fort.t' + self.frameString(frame_number)
        return open(time_file_name, 'r')
        
    
    def cell_centers(self):

        x = numpy.linspace(self.x_low+self.dx/2, self.x_high-self.dx/2, self.mx)
        y = numpy.linspace(self.y_low+self.dy/2, self.y_high-self.dy/2, self.my)

        I = numpy.indices( (self.mx, self.my) )
        
        return x[I[0]], y[I[1]]


    def cell_vertices(self):
        x = numpy.linspace(self.x_low, self.x_high, self.mx+1)
        y = numpy.linspace(self.y_low, self.y_high, self.my+1)

        I = numpy.indices( (self.mx+1, self.my+1) )
                
        return x[I[0]], y[I[1]]
        
        
    def SetFrame(self, frame_number):
    
        try:
            data_file = self.dataFile(frame_number)
            time_file = self.timeFile(frame_number)
        except:
            print "Error: could not open data or time file for frame number " + str(frame_number)
        else:
            self.frame_number = frame_number
            self.t = float( time_file.readline().split()[0] )
            time_file.close()
    
            for i in range(9):
                data_file.readline()
    
            for iy in range(self.my):

                for ix in range(self.mx):
                    self.q[ix,iy,:] = numpy.array( map( float, data_file.readline().split() ) )

                data_file.readline()
                
            data_file.close()
                
                
    def stepFrame(self):
        
        self.setFrame(self.frame_number+1)
