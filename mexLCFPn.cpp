#include "mex.h"
#include <string.h>
#include <math.h>
#include <stack>

int Cn[16] = {0};
const double pi = 3.141592653589793;

typedef unsigned int uint; //这里的int可能是64位的，MATLAB里说mxUINT32对应这个
typedef struct _Point
{
	uint x : 16; //16bit应该足够保存一般图像的下标
	uint y : 16;
} Point;

class Stack
{
  private: //member
	Point *data;
	uint index;
	uint size;
	const uint sz = sizeof(Point);

  public:
	Stack(uint size = 512)
	{
		this->size = size;
		index = 0;
		data = (Point *)mxCalloc(size, sz);
	}
	~Stack() { mxFree(data); }
	inline void push(uint x, uint y)
	{
		data[index].x = x;
		data[index].y = y;
		index++;
		if (index >= size)
			expand();
	}
	inline void pop(uint &x, uint &y)
	{
		if (index == 0)
		{
			mexErrMsgTxt("Stack is empty!!!");
			return;
		}
		index--;
		x = data[index].x;
		y = data[index].y;
	}
	inline bool isempty() { return this->index == 0; }

  private:
	inline void expand()
	{
		void *p = mxCalloc(size * 2, sz);
		memcpy(p, data, size * sz);
		mxFree(data);
		data = (Point *)p;
		size *= 2;
	}
};

class Plist
{
  private: //member
	Point *data;
	uint head;
	uint tail;
	uint size;
	const uint sz = sizeof(Point);

  public:
	Plist(uint size = 1024)
	{
		this->size = size;
		head = 0;
		tail = 0;
		data = (Point *)mxCalloc(size, sz);
	}
	~Plist() { mxFree(data); }
	inline void reset()
	{
		head = 0;
		tail = 0;
	}
	inline void push(uint x, uint y)
	{
		data[tail].x = x;
		data[tail].y = y;
		tail++;
		if (tail >= size)
			expand();
	}
	inline bool isgrowing() { return head < tail; }
	inline void grow() { head++; }
	inline uint getx() { return data[head].x; }
	inline uint gety() { return data[head].y; }
	inline void start() { head = 0; }
	inline bool notend() { return head < tail; }

  private:
	inline void expand()
	{
		void *p = mxCalloc(size * 2, sz);
		memcpy(p, data, size * sz);
		mxFree(data);
		data = (Point *)p;
		size *= 2;
	}
};

#define TRUE 4
#define FALSE 0
#define SKIP 1
#define MASK 2

const double th = 1e-5;

// floodfill时查找邻接的顺序 {row, col}
int nb[8][2] = {{-1, 0}, {-1, -1}, {0, -1}, {1, -1}, {1, 0}, {1, 1}, {0, 1}, {-1, 1}};

void suppress(double *I, char *M, uint rows, uint cols)
{
	//MATLAB里内存是列向为第一维
	uint r, c, s, x, rows_1 = rows - 1, rows_2 = rows - 2, counter = 0;
	uint pos = 0;
	uint isolated = 0;
	Stack stk;

	double *Ib = I;			//前一列列首
	double *Ic = I + rows;  //当前列列首
	double *In = Ic + rows; //下一列列首
	char *Mc = M + rows;
	char *Mn = Mc + rows;
	int ismax = 1;

	for (c = 1; c < cols - 1; ++c)
	{
		// 这里是按行查找
		// r 表示行的下标
		r = 1;
		s = rows_1;
		while (r < rows_2)
		{
			if (Mc[r] & SKIP)
			{
				r++;
				continue;
			}
			if (Ic[r] <= Ic[r + 1])
			{
				if (Ic[r] != Ic[r + 1])
					s = r + 1;
				else if (Ic[r] > Ic[r - 1])
					s = r;
				while (r++ < rows_2 && Ic[r] <= Ic[r + 1])
					if (Ic[r] != Ic[r + 1])
						s = r + 1;
			}
			else
			{
				if (Ic[r] <= Ic[r - 1])
				{
					r += 2;
					continue;
				}
				s = r;
			}
			if (s < rows_1)
			{
				ismax = 1;
				counter = 0;
				isolated = 0;
				for (x = s - 1; x <= r + 1; ++x)
				{
					if (Ic[r] >= In[x])
					{
						if (Ic[r] == In[x])
							counter++;
						else
						{
							Mn[x] |= SKIP;
							if (In[x] < th)
								isolated++;
						}
					}
					else
						ismax = 0;
					if (Ic[r] >= Ib[x])
					{
						if (Ic[r] == Ib[x])
							counter++;
						else if (Ib[x] < th)
							isolated++;
					}
					else
						ismax = 0;
				}
				if (Ic[r - 1] < th)
					isolated++;
				if (Ic[r + 1] < th)
					isolated++;
				if (ismax && isolated < 6)
				{
					if (counter > 0 || s < r)
					{
						stk.push(c, r);
					}
					while (s <= r)
						Mc[s++] |= TRUE;
				}
				s = rows;
			}
			r += 2;
		}
		Mc[rows_1] = FALSE;
		Ib += rows;
		Ic += rows;
		In += rows;
		Mc += rows;
		Mn += rows;
	}
	// 除去假的局部极值

	// 用floodfill找每个局部极大值周围是否有相同的值
	Plist plist(512);
	uint tr, tc, ind = 0;
	double k = 0;
	char u;
	while (!stk.isempty())
	{
		stk.pop(c, r);
		if (M[c * rows + r] & MASK)
			continue;
		plist.reset();
		plist.push(c, r);
		ismax = 1;
		while (plist.isgrowing())
		{
			r = plist.gety();
			c = plist.getx();
			k = I[c * rows + r];
			for (s = 0; s < 8; s++)
			{
				tr = r + nb[s][0];
				tc = c + nb[s][1];
				if (tr >= rows || tc >= cols)
					continue;
				ind = tc * rows + tr;
				if (!(M[ind] & MASK) && I[ind] == k)
				{
					plist.push(tc, tr);
					if (!(M[ind] & TRUE))
						ismax = 0;
					M[ind] |= MASK;
				}
			}
			plist.grow();
		}
		u = ismax == 0 ? MASK : TRUE | MASK;
		for (plist.start(); plist.notend(); plist.grow())
		{
			M[plist.getx() * rows + plist.gety()] = u;
		}
	}
}

void LCFP(double *S, double *O, double sigma, int H, int W, int radius)
{
	int x, y, offset, k, xH;
	double cs[16] = {0}; //采样缓存
	double A[4], B[8], f1, f2, R;
	double a, b, c, d;

	double cosp8 = cos(pi / 8);
	double sinp8 = sin(pi / 8);
	double sqrt2 = sqrt(2) / 2; //乘以 sqrt(2)/2 比除以 sqrt(2)快得多
	// p = 5*10^-7
	double threshold = 2 * 5 * log(10) * 8 * sigma * sigma;
	// double threshold = 10*sqrt(2)*sigma;
	// 0 1 2 3 4 5 6 7 8
	// 圆环上点的偏移量都是固定的，这里先计算，起到加速作用
	// 注意 MATLAB 是列先的， C是行先的。 这里应该用列先

	// x: col, y: row
	for (x = radius; x < W - radius; x++)
	{
		xH = x * H;
		for (y = radius; y < H - radius; y++)
		{
			offset = xH + y;

			// 圆环采样
			for (k = 0; k != 16; k++)
				cs[k] = S[offset + Cn[k]];

			// 为了加速计算FFT16，这里用了一些缓存方法
			// A = cs[1:4]+cs[9:12]-cs[5:8]-cs[13:16]
			// B = cs[1:8]-cs[9:16]
			for (k = 0; k < 4; k++)
			{
				a = cs[k];
				b = cs[k + 4];
				c = cs[k + 8];
				d = cs[k + 12];
				A[k] = a + c - b - d;
				B[k] = a - c;
				B[k + 4] = b - d;
			}
			// a: f1c, b: f1s, c: f2c, d: f2s
			a = B[0] + cosp8 * (B[1] - B[7]) + sinp8 * (B[3] - B[5]) + (B[2] - B[6]) * sqrt2;
			b = -B[4] - sinp8 * (B[1] + B[7]) - cosp8 * (B[3] + B[5]) - (B[2] + B[6]) * sqrt2;
			c = A[0] + (A[1] - A[3]) * sqrt2;
			d = -A[2] - (A[1] + A[3]) * sqrt2;
			f1 = a * a + b * b;
			f2 = c * c + d * d;
			R = f2 - f1;
			if (R > threshold)
				O[offset] = R;
		}
	}
}

void calc_offset(int R, int H)
{
	int x, y;
	double t = 0.0;
	double dt = 2*pi/16.0;
	for (int i = 0; i < 16; i++)
	{
		x = round(R * cos(t));
		y = round(R * sin(t));
		Cn[i] = x * H + y;
		//printf("x:%2d, y:%2d\n", x, y);
		t += dt;
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 3)
		mexErrMsgTxt("Three parameters were needed.");
	if (!(mxIsDouble(prhs[0])))
		mexErrMsgTxt("Function only support double type");
	if (nlhs != 2)
		mexErrMsgTxt("Function will output both response and local maximum");
	mwSize nd = mxGetNumberOfDimensions(prhs[0]);
	if (nd != 2)
		mexErrMsgTxt("Function only support 2d matrix");
	const mwSize *sz = mxGetDimensions(prhs[0]);
	// create response map
	plhs[0] = mxCreateDoubleMatrix(sz[0], sz[1], mxREAL);
	// create a temp map to assist nms
	plhs[1] = mxCreateNumericMatrix(sz[0], sz[1], mxINT8_CLASS, mxREAL);
	// get matrix address
	double *out = (double *)mxGetPr(plhs[0]);
	char *mask = (char *)mxGetData(plhs[1]);
	double *src = (double *)mxGetPr(prhs[0]);
	// get standard deviation
	double sigma = *((double *)mxGetPr(prhs[1]));
	double radius = *((double *)mxGetPr(prhs[2]));
	calc_offset((int)radius, sz[0]);
	LCFP(src, out, sigma, sz[0], sz[1], (int)radius);
	suppress(out, mask, sz[0], sz[1]);
	//mxFree(sz);
}
