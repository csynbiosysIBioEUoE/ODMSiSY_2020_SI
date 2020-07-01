#By Tim Sheerman-Chase 2016
#Released under the CC0 license
from __future__ import print_function
import itertools, copy

class RunningAverage(object):
	def __init__(self):
		self.Val = 0.0
		self.Count = 0

	def Update(self, valIn):
		scaling = 1.0 / float(self.Count + 1)
		self.Val = float(valIn) * scaling + self.Val * (1.0 - scaling)
		self.Count+=1
		return self.Val

	def Get(self):
		return self.Val

def TestRunningAverage():
	ra = RunningAverage()
	answers = [0.5,17.75,12.0666666667,7.55,6.84,6.7666666667,6.5142857143,6.65,6.1888888889,
		5.69,4.6272727273,4.4083333333,4.0692307692,3.7785714286,3.66,3.68125,3.8176470588]
	for val, answer in zip([0.5,35,0.7,-6,4,6.4,5,7.6,2.5,1.2,-6,2,0,0,2,4,6], answers):
		result = ra.Update(val)
		if abs(result - answer) > 1e-6:
			raise RuntimeError("Error in running average")		

def SampStdDev(vals, valsSum = None):
	if valsSum is None:
		valsSum = sum(vals)
	me = float(valsSum) / len(vals)
	sq = [(v - me)**2.0 for v in vals]	
	sqSum = sum(sq)
	return (sqSum / (len(vals)-1))**0.5
	
def SampPearsonCorrel(a, b, asum=None, bsum=None, adev=None, bdev=None):
	if len(a) != len(b):
		raise ValueError("Input lengths must match")
	if asum is None:
		asum = sum(a)
	if bsum is None:
		bsum = sum(b)
	if adev is None:
		adev = SampStdDev(a, asum)
	if bdev is None:
		bdev = SampStdDev(b, bsum)
	
	terms = [s1*s2 for s1, s2 in zip(a, b)]
	return (sum(terms) - asum * bsum / len(a)) / ((len(a)-1) * adev * bdev)

def isiterable(data):
	try:
		chkit = iter(data)
	except TypeError, te:
		return False
	return True

def array(data):
	arr = SimpleArray()
	arr.shape = []
	if not isiterable(data):
		raise ValueError("Input must be iterable")
	arr._ValidateShape(data, 0)
	arr.data = data
	arr.shape = tuple(arr.shape)
	return arr

def zeros(shape):
	arr = SimpleArray()
	arr.shape = tuple(shape)
	if not isiterable(shape):
		raise ValueError("Shape must be iterable")
	arr.data = []
	def GenData(shape, leafs):
		dim = shape.pop(0)
		newLeafs = []
		if len(shape) == 0:
			for l in leafs:
				for i in range(dim):
					l.append(0.0)
		else:
			for l in leafs:
				for i in range(dim):
					nl = []
					l.append(nl)
					newLeafs.append(nl)
			GenData(shape, newLeafs)
			
	GenData(list(arr.shape)[:], [arr.data])
	return arr

def identity(dim):
	mat = zeros((dim, dim))
	for i in range(dim):
		mat.data[i][i] = 1.0
	return mat

class SimpleArray(object):
	def __init__(self):
		self.shape = []
		self.data = None

	def _ValidateShape(self, data, ind):
		if not isiterable(data):
			return
		if ind >= len(self.shape):
			self.shape.append(len(data))
		else:
			if self.shape[ind] != len(data):
				raise ValueError("Inconsistent shape")
		
		for val in data:
			self._ValidateShape(val, ind+1)

	def __repr__(self):
		return (self.__class__.__name__+"("+str(self.data)+")")

	def dot(self, rhs):
		#This also standard matrix multiplication
		if not isinstance(rhs, self.__class__) and isiterable(rhs):		
			rhs = array(rhs)
		if len(self.shape) != 2 or len(rhs.shape) != 2:
			raise NotImplementedError("Only implemented for 2D matricies")
		if self.shape[1] != rhs.shape[0]:
			raise ValueError("Matrix size mismatch")
		m = self.shape[1]
		result = zeros((self.shape[0], rhs.shape[1]))
		for i in range(result.shape[0]):
			for j in range(result.shape[1]):
				tot = 0.0
				for k in range(m):
					aik = self.data[i][k]
					bkj = rhs.data[k][j]
					tot += aik * bkj
				result.data[i][j] = tot

		return result

	def __add__(self, rhs):
		if not isinstance(rhs, self.__class__) and isiterable(rhs):		
			rhs = array(rhs)
		if self.shape != rhs.shape:
			raise ValueError("Matrix size mismatch")
		result = zeros(self.shape)
		def AddFunc(ind, a, b, out):
			if ind == len(self.shape):
				raise ValueError("Invalid matrix")
			if ind == len(self.shape) - 1:
				for i, (av, bv) in enumerate(zip(a, b)):
					out[i] = av + bv
			else:
				ind += 1
				for ar, br, outr in zip(a, b, out):
					AddFunc(ind, ar, br, outr)		

		AddFunc(0, self.data, rhs.data, result.data)
		return result

	def __mul__(self, rhs):
		if not isinstance(rhs, self.__class__) and isiterable(rhs):		
			rhs = array(rhs)

		if isinstance(rhs, self.__class__):
			#Do element-wise matrix multiplication

			if self.shape != rhs.shape:
				raise ValueError("Matrix size mismatch")
			result = zeros(self.shape)
			def MultFunc(ind, a, b, out):
				if ind == len(self.shape):
					raise ValueError("Invalid matrix")
				if ind == len(self.shape) - 1:
					for i, (av, bv) in enumerate(zip(a, b)):
						out[i] = av * bv
				else:
					ind += 1
					for ar, br, outr in zip(a, b, out):
						MultFunc(ind, ar, br, outr)		

			MultFunc(0, self.data, rhs.data, result.data)
			return result
		else:
			#Do scalar multiplication to entire matrix (element-wise)
			result = zeros(self.shape)
			def MultFunc2(ind, a, b, out):
				if ind == len(self.shape):
					raise ValueError("Invalid matrix")
				if ind == len(self.shape) - 1:
					for i, av in enumerate(a):
						out[i] = av * b
				else:
					ind += 1
					for ar, outr in zip(a, out):
						MultFunc2(ind, ar, b, outr)		

			MultFunc2(0, self.data, rhs, result.data)
			return result

	def __sub__(self, rhs):
		if not isinstance(rhs, self.__class__) and isiterable(rhs):		
			rhs = array(rhs)
		if self.shape != rhs.shape:
			raise ValueError("Matrix size mismatch")
		negRhs = rhs * -1
		return self + negRhs

	def conj(self):
		result = zeros(self.shape)
		def AddFunc(ind, a, out):
			if ind == len(self.shape):
				raise ValueError("Invalid matrix")
			if ind == len(self.shape) - 1:
				for i, av in enumerate(a):
					if isinstance(av, complex):
						out[i] = av.conjugate()
					else:
						out[i] = av
			else:
				ind += 1
				for ar, outr in zip(a, out):
					AddFunc(ind, ar, outr)		

		AddFunc(0, self.data, result.data)
		return result

	@property
	def T(self):
		if len(self.shape) <= 1:
			return self
		if len(self.shape) != 2:
			raise NotImplementedError("Only implemented for 1D and 2D matricies")	
		result = zeros((self.shape[1], self.shape[0]))
		for i in range(self.shape[0]):
			for j in range(self.shape[1]):
				result.data[j][i] = self.data[i][j]
		return result

	def copy(self):
		return array(copy.deepcopy(self.data))

def perm_parity(lst):
	'''\
	Given a permutation of the digits 0..N in order as a list, 
	returns its parity (or sign): +1 for even parity; -1 for odd.
	From https://code.activestate.com/recipes/578227-generate-the-parity-or-sign-of-a-permutation/
	'''
	parity = 1
	for i in range(0,len(lst)-1):
		if lst[i] != i:
			parity *= -1
			mn = min(range(i,len(lst)), key=lst.__getitem__)
			lst[i],lst[mn] = lst[mn],lst[i]
	return parity	

def det(mat):
	if not isinstance(mat, SimpleArray) and isiterable(mat):		
		mat = array(mat)
	if len(mat.shape) != 2:
		raise NotImplementedError("Only implemented for 2D matricies")
	if mat.shape[0] != mat.shape[1]:
		raise ValueError("Matrix must be square")
	n = mat.shape[0]
	
	#Leibniz formula for the determinant
	total2 = 0.0
	for perm in itertools.permutations(range(n)):
		total1 = 1.0
		for i, j in zip(range(n), perm):
			total1 *= mat.data[i][j]
		total2 += perm_parity(list(perm)) * total1
	return total2

def delete(mat, ind, axis=0):
	if not isinstance(mat, SimpleArray) and isiterable(mat):		
		mat = array(mat)
	newShape = list(mat.shape)[:]
	newShape[axis] -= len(ind)
	result = zeros(newShape)
	def CopyWithDelete(currentAxis, ind, axis, matdata, resultData):
		if currentAxis < len(mat.shape)-1:
			count = 0
			for i, data in enumerate(matdata):
				if currentAxis == axis and i in ind:
					continue
				CopyWithDelete(currentAxis + 1, ind, axis, data, resultData[count])
				count += 1
		else:
			count = 0
			for i, val in enumerate(matdata):
				if currentAxis == axis and i in ind:
					continue
				resultData[count] = val
				count += 1

	CopyWithDelete(0, ind, axis, mat.data, result.data)
	return result

def adjugate(mat):
	if not isinstance(mat, SimpleArray) and isiterable(mat):		
		mat = array(mat)
	#Find the adjugate matrix
	if len(mat.shape) != 2:
		raise NotImplementedError("Only implemented for 2D matricies")
	result = zeros(mat.shape)
	for i in range(mat.shape[0]):
		for j in range(mat.shape[1]):
			submat = delete(mat, [j], 0)
			submat2 = delete(submat, [i], 1)
			result.data[i][j] = ((-1)**(i+j))*det(submat2)
	return result

def inv_by_adjugate(mat, eps=1e-8):
	#Find inverse based on find the adjugate matrix
	#which is inefficient for large matrices.
	if not isinstance(mat, SimpleArray) and isiterable(mat):		
		mat = array(mat)
	if len(mat.shape) != 2:
		raise NotImplementedError("Only implemented for 2D matricies")
	if mat.shape[0] != mat.shape[1]:
		raise ValueError("Matrix must be square")
	mdet = det(mat)
	if abs(mdet) < eps:
		raise RuntimeError("Matrix is not invertible (its determinant is zero)")
	return adjugate(mat) * (1.0 / float(mdet))

def inv_by_gauss_jordan(mat, eps=1e-8):
	#Find inverse based on gauss-jordan elimination.

	if not isinstance(mat, SimpleArray) and isiterable(mat):		
		mat = array(mat)
	if len(mat.shape) != 2:
		raise NotImplementedError("Only implemented for 2D matricies")
	if mat.shape[0] != mat.shape[1]:
		raise ValueError("Matrix must be square")
	mdet = det(mat)
	if abs(mdet) < eps:
		raise RuntimeError("Matrix is not invertible (its determinant is zero)")

	#Create aux matrix
	n = mat.shape[0]
	auxmat = identity(n)

	#Convert to echelon (triangular) form
	mat = copy.deepcopy(mat)
	for i in range(n):
		#Find highest value in this column
		maxv = 0.0
		maxind = None
		for r in range(i, n):
			v = mat.data[r][i]
			if maxind is None or abs(v) > maxv:
				maxv = abs(v)
				maxind = r
		
		if maxind != i:
			#Swap this to the current row, for numerical stability
			mat.data[i], mat.data[maxind] = mat.data[maxind], mat.data[i]
			auxmat.data[i], auxmat.data[maxind] = auxmat.data[maxind], auxmat.data[i]

		activeRow = mat.data[i]
		activeAuxRow = auxmat.data[i]
		for r in range(i+1, n):
			scale = float(mat.data[r][i]) / float(mat.data[i][i])
			cursorRow = mat.data[r]
			cursorAuxRow = auxmat.data[r]
			for c in range(n):
				cursorRow[c] -= scale * activeRow[c]
				cursorAuxRow[c] -= scale * activeAuxRow[c]

	#Finish elimination
	for i in range(n-1, -1, -1):
		activeRow = mat.data[i]
		activeAuxRow = auxmat.data[i]
		for r in range(i, -1, -1):
			cursorRow = mat.data[r]
			cursorAuxRow = auxmat.data[r]
			if r == i:
				scaling = activeRow[i]
				for c in range(n):
					cursorRow[c] /= scaling
					cursorAuxRow[c] /= scaling
			else:
				scaling = cursorRow[i]
				for c in range(n):
					cursorRow[c] -= activeRow[c] * scaling
					cursorAuxRow[c] -= activeAuxRow[c] * scaling
			
	return auxmat

def pinv(mat):
	#Find pseudo-inverse of a matrix
	if not isinstance(mat, SimpleArray) and isiterable(mat):		
		mat = array(mat)
	if len(mat.shape) != 2:
		raise NotImplementedError("Only implemented for 2D matricies")
	matconjh = mat.conj().T
	return inv_by_adjugate(matconjh.dot(mat)).dot(matconjh)

def TestArray():
	data = [1,2,3]
	a = array(data)
	if a.shape != (3,):
		raise RuntimeError("Incorrect shape")

	data = [[1,2,3], [4,5,6]]
	a = array(data)
	if a.shape != (2, 3):
		raise RuntimeError("Incorrect shape")

	ok = False
	try:
		a = array([[1,2,3], [4,5]])
	except ValueError:
		ok = True
	if not ok:
		raise RuntimeError("Failed to detect invalid input")

	a = array([[1,2],[3,4],[5,6]])
	b = array([[7,8,9],[10,11,12]])
	result = a.dot(b)
	ans= array([[ 27,  30,  33],
	   [ 61,  68,  75],
	   [ 95, 106, 117]])
	err = result - ans
	for row in err.data:
		for errVal in row:
			if abs(errVal) > 1e-6:
				raise RuntimeError("Matrix dot product has wrong result")

	data = array([1,2+4j,3])
	result = data.conj()
	err = result - array([1,2-4j,3])
	for errVal in err.data:
		if abs(errVal) > 1e-6:
			raise RuntimeError("Matrix conjugation has wrong result")

	data = array([[7,2,3], [12,5,6], [7,8,9]])
	result = det(data)
	if abs(result - 30.0) > 1e-6:
			raise RuntimeError("Matrix determinant has wrong result")

	data = array([[12,1,2,3], [5,4,5,6], [11,7,8,9], [5,6,2,5]])
	result = det(data)
	if abs(result + 273.0) > 1e-6:
			raise RuntimeError("Matrix determinant has wrong result")
	
	data = array([[7,2,3], [12,5,6], [7,8,9]])
	ans = array([[7, 2, 3], [7, 8, 9]])
	result = delete(data, [1], 0)
	err = result - ans
	for row in err.data:
		for errVal in row:
			if abs(errVal) > 1e-6:
				raise RuntimeError("Matrix delete has wrong result")
	
	data = array([[7,2,3], [12,5,6], [7,8,9]])
	ans = array([[2],[5],[8]])
	result = delete(data, [0, 2], 1)
	err = result - ans
	for row in err.data:
		for errVal in row:
			if abs(errVal) > 1e-6:
				raise RuntimeError("Matrix delete has wrong result")

	data = array([[7,2,3], [12,5,6], [7,8,9]])
	result = adjugate(data)
	ans = array([[-3,6,-3],[-66,42,-6],[ 61,-42,11]])
	err = result - ans
	for row in err.data:
		for errVal in row:
			if errVal > 1e-6:
				raise RuntimeError("Matrix adjugate has wrong result")

	data = array([[12,1,2,3], [5,4,5,6], [11,7,8,9], [5,6,2,5]])
	result = adjugate(data)
	ans = array([[-21,42,-21,0],[48,190,-121,-39],[51,145,-157,78],[-57,-328,229,-39]])
	err = result - ans
	for row in err.data:
		for errVal in row:
			if abs(errVal) > 1e-6:
				raise RuntimeError("Matrix adjugate has wrong result")

	data = array([[7,2,3], [12,5,6], [7,8,9]])
	result = inv_by_adjugate(data)
	ans = array([[-0.1,0.2,-0.1],[-2.2,1.4,-0.2],[2.03333333,-1.4,0.36666667]])
	err = result - ans
	for row in err.data:
		for errVal in row:
			if abs(errVal) > 1e-6:
				raise RuntimeError("Matrix inverse has wrong result")

	data = array([[12,1,2,3], [5,4,5,6], [11,7,8,9], [5,6,2,5]])
	result = inv_by_adjugate(data)
	ans = ([[ 0.07692308, -0.15384615,  0.07692308,  0.		],
		[-0.17582418, -0.6959707 ,  0.44322344,  0.14285714],
		[-0.18681319, -0.53113553,  0.57509158, -0.28571429],
		[ 0.20879121,  1.2014652 , -0.83882784,  0.14285714]])
	err = result - ans
	for row in err.data:
		for errVal in row:
			if abs(errVal) > 1e-6:
				raise RuntimeError("Matrix inverse has wrong result")

	ok = False
	data = array([[12,1,2,3], [0,0,0,0], [11,7,8,9], [5,6,2,5]])
	try:
		result = inv_by_adjugate(data)
	except RuntimeError:
		ok = True
	if not ok:
		raise RuntimeError("Failed to detect uninverible matrix")
	
	data = array([[7,2], [12,5], [7,8]])
	result = pinv(data)
	ans = [[ 0.0697467 ,  0.08312522, -0.06938994],
		[-0.07599001, -0.06243311,  0.18301819]]
	err = result - ans
	for row in err.data:
		for errVal in row:
			if abs(errVal) > 1e-6:
				raise RuntimeError("Matrix pseudo-inverse has wrong result")

	data = [[2,-1,0],[-1,2,-1],[0,-1,2]]
	result = inv_by_gauss_jordan(data)
	ans = array([[0.75, 0.5, 0.25], [0.5, 1.0, 0.5], [0.25, 0.5, 0.75]])
	err = result - ans
	for row in err.data:
		for errVal in row:
			if abs(errVal) > 1e-6:
				raise RuntimeError("Matrix inverse has wrong result")

	data = array([[7,2,3], [12,5,6], [7,8,9]])
	result = inv_by_gauss_jordan(data)
	ans = array([[-0.1,0.2,-0.1],[-2.2,1.4,-0.2],[2.03333333,-1.4,0.36666667]])
	err = result - ans
	for row in err.data:
		for errVal in row:
			if abs(errVal) > 1e-6:
				raise RuntimeError("Matrix inverse has wrong result")

	data = array([[12,1,2,3], [5,4,5,6], [11,7,8,9], [5,6,2,5]])
	result = inv_by_gauss_jordan(data)
	ans = ([[ 0.07692308, -0.15384615,  0.07692308,  0.		],
		[-0.17582418, -0.6959707 ,  0.44322344,  0.14285714],
		[-0.18681319, -0.53113553,  0.57509158, -0.28571429],
		[ 0.20879121,  1.2014652 , -0.83882784,  0.14285714]])
	err = result - ans
	for row in err.data:
		for errVal in row:
			if abs(errVal) > 1e-6:
				raise RuntimeError("Matrix inverse has wrong result")

if __name__=="__main__":
	TestRunningAverage()
	TestArray()

	print ("All tests passed")
