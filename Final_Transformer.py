#Name: Eleftherios Vganges
#Date: 09/29/2023
#Program: The purpose of this program is to modify the ray tracer to include global and local coordinate systems that group objects under parents.
#           This allows for translations and rotations to be done on the object which reposition its coordinates.
''' credits
[1] https://observablehq.com/@infowantstobeseen/basic-ray-casting  
[2] http://lousodrome.net/blog/light/2020/07/03/intersection-of-a-ray-and-a-plane/
[3] http://www.bmsc.washington.edu/people/merritt/graphics/quadrics.html
[4] Fundamentals of Computer Graphics 5th ed, Steve Marschner Peter Shirley
for assignment 3
[5] https://observablehq.com/@cse4413msstate/basic-transformations
[6] https://observablehq.com/@cse4413msstate/understanding-transform-composition
'''
import math
import numpy as np
from PIL import Image, ImageOps

def matrixer(object, matrix):   # main matrix multiplication function, accepts object being moved and matrix of movement
    # for sphere take center point and relocate based on output of matrix multiplication
    if object.id == 1:
        spheremat = np.matrix([[object.center.x], [object.center.y], [object.center.z], [1]])
        out = np.matmul(matrix, spheremat)
        object.center = point(out.item(0,0), out.item(1,0), out.item(2,0))

    # for triangle, does the same but to each vertex
    if object.id == 2:
        vertextrans = np.matrix([[object.a.x], [object.a.y], [object.a.z], [1]])
        out = np.matmul(matrix, vertextrans)
        object.a = point(out.item(0,0), out.item(1,0), out.item(2,0))
        
        vertextrans = np.matrix([[object.b.x], [object.b.y], [object.b.z], [1]])
        out = np.matmul(matrix, vertextrans)
        object.b = point(out.item(0,0), out.item(1,0), out.item(2,0))

        vertextrans = np.matrix([[object.c.x], [object.c.y], [object.c.z], [1]])
        out = np.matmul(matrix, vertextrans)
        object.c = point(out.item(0,0), out.item(1,0), out.item(2,0))

    # for the coordinate systems, allows their position to be updated when moved.
    if object.id == 0:
        center = object.mat
        out = np.matmul(matrix, center)
        object.update(out)

    # for the camera, recenters eye based on output of multiplication
    if object.id == 4:
        eye = np.matrix([[object.e.x],[object.e.y],[object.e.z],[1]])
        out = np.matmul(matrix, eye)
        object.e = point(out.item(0,0), out.item(1,0), out.item(2,0))


def RecursiveAlgorithm(subject, objlist, matrices): # recursive algorithm that sends an object through the multiplication process, and affects each child 
    matrixer(subject, matrices) # moves coordinate system
    for j in objlist:
        if j.parent == subject:
            matrixer(j, matrices)   # moves object

    if subject.child != None:   # if object has children, calls function again on each child
        for i in subject.child:
            RecursiveAlgorithm(i, objlist, matrices)
        
    

def main():

    fname = input("Filename: ")  # input of file to read from
    
    op = 0
    while op != int(1) and op != int(2):    # choice of perspective or orthographic
        op = int(input("Orthographic-1 or Perspective-2: "))
    
    itemlist = reader(fname)    # reads items to a list and assigns them below in (width, height), image plane, cam, and the list of objects
    
    width = itemlist[0][0]
    height = itemlist[0][1]
    imp = itemlist[1]
    cam = itemlist[2]
    objlist = itemlist[3]

    systems = itemlist[4] # each coordinate system
    transforms = itemlist[5]    # each transform to be done
    
    image = Image.new("RGB", (width, height), (255, 255, 255))  # creates a new image with width and height and loads pixels

    imnum = 0
    for x in transforms:    # for each transformation in the list, output an image
        pixels = image.load()
        
        transobj = x.subject    # finds object to be moved and either sends it directly if camera
        if isinstance(transobj, Camera):
            matrixer(transobj, x.matrices)
        else:   # or sends it through recursion if an coordinate system
            RecursiveAlgorithm(transobj, objlist, x.matrices)
        print("Transforming by", x)

        for x in range(width):
            for y in range(height): # for each pixel

                u = imp.l + (imp.r - imp.l) * (x + 0.5) / width     # calculates u and v aspect ratio for image projection (4.1) from textbook
                v = imp.b + (imp.t - imp.b) * (y + 0.5) / height

                if op == 1:     # orthographic proj using formula on page 85
                    shotray = ray(cam.e + cam.u.mult(u) + cam.v.mult(v), cam.w.mult(-1))

                elif op == 2:   # perspective proj using formula on page 86
                    dist = np.sqrt((u - cam.e.x) ** 2 + (v - cam.e.y) ** 2 + (imp.z - cam.e.z) ** 2) # dist = |o(center of pixel) - e|    o = (u , v, imp.z)
                    shotray = ray(cam.e, cam.w.mult(dist * -1) + cam.u.mult(u) + cam.v.mult(v)) 

                prvdata = [False, np.inf]   # sets a holder through iterations
                for z in objlist:   # for each object
                    hit_data = z.intersect(shotray) #gets hit data, [T/F, int] 

                    if hit_data[0] == True: #if it does intersect

                        if hit_data[1] < prvdata[1] and hit_data[1] >= 0: #if the shapes intersection is less than any other the other shapes
                            prvdata = hit_data
                            pixels[x,y] = z.color   #colors pixel of smallest shape

    
        image = ImageOps.flip(image)    # flip because of y being defined downwards 
        imnum = imnum + 1
        image.save(str(imnum) + "OutputImage.png")   
        image.show()

class Translate(object): # creates a translation matrix with a given x y and z
    def __init__(self, x, y, z):
        self.mat = np.matrix([[1,0,0,x], [0,1,0,y], [0,0,1,-1 * z], [0,0,0,1]])

class Rotate(object): # creates a rotation matrix about either x y or z and bases it on given theta 
    def __init__(self, type, theta):
        if type == 0: #x
            self.mat = np.matrix([[1,0,0,0], [0, math.cos(theta), (-1 * math.sin(theta)), 0], [0, math.sin(theta), math.cos(theta), 0], [0,0,0,1]])
        elif type == 1: #y
            self.mat = np.matrix([[math.cos(theta), 0, math.sin(theta), 0], [0,1,0,0], [(-1 * math.sin(theta)), 0, math.cos(theta), 0], [0,0,0,1]])
        elif type == 2: #z
            self.mat = np.matrix([[math.cos(theta),(-1 * math.sin(theta)),0,0], [math.sin(theta), math.cos(theta), 0, 0], [0,0, 1, 0], [0,0,0,1]])
        else:
            print("Incorrect Type\n")
            exit()


class point(object):    # creates a simple object that holds 3 values
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def copy(self):
        return self

    #supports addition subtraction dot product multiplication and subtraction by a value normalization and cross product
    def __add__(self, other):
        return point(self.x + other.x, self.y + other.y, self.z + other.z)
    
    def __sub__(self, other):
        return point(self.x - other.x, self.y - other.y, self.z - other.z)
    
    def dot(p1, p2):
        return p1.x * p2.x + p1.y * p2.y + p1.z * p2.z
    
    def mult(self, int):
        return point(self.x * int, self.y * int, self.z * int)
    
    def floatsub(p1, float):
        return point(p1.x - float, p1.y - float, p1.z - float)
    
    def normalize(self):
        mg = np.sqrt(self.x * self.x, self.y * self.y, self.z * self.z)
        self.x = self.x / mg
        self.y = self.y / mg
        self.z = self.z / mg

    def cross(self, p1):
        x = [self.x, self.y, self.z]
        y = [p1.x, p1.y, p1.z]
        return np.cross(x, y)

class CoordSys(object): # coordsystem object which accepts 3 vectors and an origin
    def __init__(self, u, v, w, origin, name = None, parent = None, child = []):
   
        self.origin = origin
        self.parent = parent
        self.child = child
        self.name = name
        # converts vectors and origin to matrix for multiplication
        self.mat = np.matrix([[u.x, u.y, u.z, origin.x], [v.x, v.y, v.z, origin.y], [w.x, w.y, w.z, origin.z], [0, 0, 0 , 1]])
        self.id = 0

    # when moved, updates matrix
    def update(self, mat):
        self.mat = mat

class Camera(object): # defines a camera with a given eye, and defaults of u,v,w, and an up
    def __init__(self, e, w = point(0,0,-1), u = point(1,0,0), v = point(0,1,0), up = point(0.00001,1,0.00001)):
        self.e = e
        self.u = u
        self.v = v
        self.w = w
        self.up = up
        self.id = 4
        #adds a camera matrix based off of the 3 vectors and its eye
        self.mat = np.matrix([[u.x, u.y, u.z, e.x], [v.x, v.y, v.z, e.y], [w.x, w.y, w.z, e.z], [0, 0, 0 , 1]])

class Transform(object): # an object that holds a transformation matrix and the object being transformed
    def __init__(self, subject, matrices):
        self.subject = subject
        self.matrices = matrices

class ImagePlane(object):   # defines an image plane from two points with leftmost rightmost topmost and bottommost, and sets the z to first point given
    def __init__(self, p1, p2):
        self.l = p1.x
        self.r = p2.x
        self.b = p1.y
        self.t = p2.y
        self.z = p1.z

class ray(object):  # an object that holds an origin point and a direction vector
    def __init__(self, origin = point(0,0,0), direction = point(0,0,0)):
        
        self.origin = origin
        self.direction = direction


class Surface(object):  # a default parent class for all objects
    def intersect(self, ray):
        return [False, -1]

class sphere(Surface):  # class for a sphere object takes in a center point radius value and a color value
    def __init__(self, center, radius, color, parent = None):
        self.center = center
        self.radius = radius
        self.color = color
        self.parent = parent
        self.id = 1

    def intersect(self, ray):
        
        distance = ray.origin - self.center # gets the rays distance from center
        
        a = point.dot(ray.direction, ray.direction) # gets basic a b and c value described in the textbook pg 87
        b = point.dot(ray.direction, distance)
        c = point.dot((distance), (distance)) - self.radius * self.radius

        discrimant = (b * b) - (a * c)  # checks to see if discriminant is 0 because if it is there is no intersection
        if discrimant < 0:
            return [False, -1]
        
        int1 = (-b - np.sqrt(b*b - a*c)) / a    # solves for the quadratic equation with the values calculated above
        int2 = (-b + np.sqrt(b*b - a*c)) / a
        
        
        if int2 < 0: # if largest point is less than 0 then the shape is behind the view plane 
            return [False, -1]
        
        elif int1 < 0:  # if largest point is greater than 0, and smallest point isnt return largest
            return [True, int2]
        
        elif int2 != int1: # if they are both hit, return the lesser value
            inthold = 0
            if int1 < int2:
                inthold = int1
            else:
                inthold = int2
            return [True, inthold]
        
        else:   # if they are equal just return first
            return [True, int1]
    
class triangle(Surface):    # defines a triangle with 3 points and a color value
    def __init__(self, p1, p2, p3, color, parent = None):
        self.a = p1
        self.b = p2
        self.c = p3
        self.color = color
        self.parent = parent
        self.id = 2


    def intersect(self, ray):   # finds the intersection point between a triangle and a ray
        # this is the matrix multiplication from section 4.4.2 wrote out in long form to save on computations
        eihf = (self.a.y - self.c.y) * ray.direction.z - (self.a.z - self.c.z) * ray.direction.y
        gfdi = (self.a.z - self.c.z) * ray.direction.x - (self.a.x - self.c.x) * ray.direction.z
        dheg = (self.a.x - self.c.x) * ray.direction.y - (self.a.y - self.c.y) * ray.direction.x
        akjb = (self.a.x - self.b.x) * (self.a.y - ray.origin.y) - (self.a.y - self.b.y) * (self.a.x - ray.origin.x)
        jcal = (self.a.z - self.b.z) * (self.a.x - ray.origin.x) - (self.a.x - self.b.x) * (self.a.z - ray.origin.z)
        blkc = (self.a.y - self.b.y) * (self.a.z - ray.origin.z) - (self.a.z - self.b.z) * (self.a.y - ray.origin.y)

        # the divisor for each equation of t, gamma, and beta, if its 0, intersection parallel and returns false to save from diving by 0
        m = (self.a.x - self.b.x) * eihf + (self.a.y - self.b.y) * gfdi + (self.a.z - self.b.z) * dheg
        if m == 0:
            return[False, -1]

        # the intersection with the triangle, if its equal to 0, return false because no intersection
        t = (-1) * ((self.a.z - self.c.z) * akjb + (self.a.y - self.c.y) * jcal + (self.a.x - self.c.x) * blkc) / m
        if t < 0:
            return [False, -1]
        
        # finds gamma, which if its not in between 0 and 1 returns false
        gamma = ((ray.direction.z * akjb) + (ray.direction.y * jcal) + (ray.direction.x * blkc)) / m #source of weird triangle mess up
        if gamma < 0 or gamma > 1:
            return [False, -1]
        
        # finds beta, which if is less than 0 or greater than 1 - gamma, returns false because intersection falls outsode
        beta = ((self.a.x - ray.origin.x) * eihf + (self.a.y - ray.origin.y) * gfdi + (self.a.z - ray.origin.z) * dheg) / m
        if beta < 0 or beta > (1 - gamma):
            return [False, -1]

        return [True, t]   # if no fail checks are true, return intersect with t as distance. again, this is the method from section 4.4.2
    
    def printtri(self):
        print("(",self.a.x,self.a.y,self.a.z,")", "(", self.b.x,self.b.y,self.b.z,")", "(", self.c.x,self.c.y,self.c.z,")", "\n")

class plane(Surface): # defines a plane as a coordinate and a normal value, source [2] used
    def __init__(self, p1, normal, color, parent = None):
        self.p1 = p1
        self.normal = normal
        self.color = color
        self.parent = parent
        self.id = 3
    
    def intersect(self, ray):   # finds intersect between point and ray
        dn = point.dot(ray.direction, self.normal)  # if 0 ray is parallel to plane
        if dn == 0:
            return [False, -1]
        
        t = point.dot((self.p1 - ray.origin), self.normal) / point.dot(ray.direction, self.normal) # finds intersection by distance * norm / direction * normal
        if t > 0: # if t is in front of image plane, return true and t
            return [True, t]
        
        else:
            return [False, -1]

def reader(filename): # custom object loader function to focus more on shape and geometry aspects
    f = open(filename, "r") # opens specified text file and reads all lines into a list
    list = f.readlines()
    finallist = []  # defines two holder lists
    objlist = []
    
    coordsyslist = [] # defines a holder list for the coord systems and transforms
    transformlist = []
    inlist = 0
    
    for x in list:  # strips all newlines from list and divides them along spaces
        x = x.rstrip('\n')
        x = x.split(" ")

        if inlist == 1: # if in the hierarchy list
            if x[0] == "h_end":
                inlist = 0
                continue
            if x[0] == "global":    # sets the required global to its given points
                coordsys = CoordSys(point(float(x[1]), float(x[2]), float(x[3])), point(float(x[4]), float(x[5]), float(x[6])), 
                           point(float(x[7]), float(x[8]), float(x[9])), point(float(x[10]), float(x[11]), float(x[12])), x[0])
                coordsyslist.append(coordsys)
            else: # for anything following the global, creates coord system and its parent
                for y in coordsyslist:
                    if y.name == x[13]:
                        coordsys = CoordSys(point(float(x[1]), float(x[2]), float(x[3])), point(float(x[4]), float(x[5]), float(x[6])), 
                           point(float(x[7]), float(x[8]), float(x[9])), point(float(x[10]), float(x[11]), float(x[12])), x[0], y)
                        coordsyslist.append(coordsys)
                        break
                    else:
                        pass

        if inlist == 2: # if in the transformation list
            if x[0] == "t_end":
                inlist = 0
                continue
            if x[1] == "translation":   # if translation, get a translation object and give it and the subject to a Transform for storage
                tmpmat = Translate(float(x[2]), float(x[3]), float(x[4]))
                if x[0] == "camera":    # if camera is moved, don't go through coord list
                    tmptransform = Transform(finallist[2], tmpmat.mat)
                    transformlist.append(tmptransform)
                else:
                    for k in coordsyslist:
                        if k.name == x[0]:
                            tmptransform = Transform(k, tmpmat.mat)
                            transformlist.append(tmptransform)
                            break
            if x[1] == "rotation":  # if rotation, creates a rotation object and gives it and the subject to a transform
                if x[2] == "x":
                    tmpmat = Rotate(0, float(x[3]))
                elif x[2] == "y":
                    tmpmat = Rotate(1, float(x[3]))
                elif x[2] == "z":
                    tmpmat = Rotate(2, float(x[3]))
                if x[0] == "camera": # if rotating camera, ignore coordsys list
                    tmptransform = Transform(finallist[2], tmpmat.mat)
                    transformlist.append(tmptransform)
                else:    
                    for k in coordsyslist:
                        if k.name == x[0]:
                            tmptransform = Transform(k, tmpmat.mat)
                            transformlist.append(tmptransform)
                            break



        if x[0] == "image": # image dimensions, set both values after into a list and add to final
            dims = [int(x[1]), int(x[2])]
            finallist.append(dims)
        if x[0] == "viewport":  # if viewport, create an imageplane object from 6 values and adds
            im = ImagePlane(point(float(x[1]), float(x[2]), float(x[3])), point(float(x[4]), float(x[5]), float(x[6])))
            finallist.append(im)
        if x[0] == "eye" or x[0] == "direction":   # if eye, create a camera with specified eye point for perspec, or sets as direction point if orthogonal
            if x[0] == "direction":
                cam = Camera(point(float(0), float(0), float(0)), point(float(x[1]), float(x[2]), float(x[3])))
            elif x[0] == "eye":
                cam = Camera(point(float(x[1]), float(x[2]), float(x[3])))
            finallist.append(cam)
        
        if x[0] == "triangle":  # if object is triangle, creates a triangle point from 9 values and 3 color values and adds to objlist
            for p in coordsyslist:
                if x[13] == p.name:
                    tri = triangle(point(float(x[1]), float(x[2]), float(x[3])), point(float(x[4]), float(x[5]), float(x[6])), 
                           point(float(x[7]), float(x[8]), float(x[9])), (int(x[10]), int(x[11]), int(x[12])), p)
                    break
            objlist.append(tri)
        
        if x[0] == "sphere":    # if object is sphere, creates sphere from 3 values for a point 1 for a radius and 3 for color and adds to objlist
            for p in coordsyslist:
                if x[8] == p.name:
                    sph = sphere(point(float(x[1]), float(x[2]), float(x[3])), float(x[4]), (int(x[5]), int(x[6]), int(x[7])), p)
                    break
            objlist.append(sph)
        
        if x[0] == "plane":     # if object is a plane, creates a plan object and sets 3 values for point 3 for normal and 3 for color
            for p in coordsyslist:
                if x[10] == p.name:
                    pln = plane(point(float(x[1]), float(x[2]), float(x[3])), point(float(x[4]), float(x[5]), float(x[6])), 
                        (int(x[7]), int(x[8]), int(x[9])), p)
                    break
            objlist.append(pln)
        
        if x[0] == "hierarchy": # enters hierarchy list
            inlist = 1
            continue
        
        if x[0] == "transforms": # enters transform list
            inlist = 2
            continue
        

    for m in coordsyslist:  # for each parent in a coord sys, sets each system labeling it as a parent to its child
        tmparr = []
        for n in coordsyslist:
            if n.parent == m:
                tmparr.append(n)
            else:
                pass
        m.child = tmparr



    finallist.append(objlist)   #adds objlist as final element of finallist and returns it
    finallist.append(coordsyslist)  # appends both lists to the final returned list
    finallist.append(transformlist)
    return(finallist)
         
main()