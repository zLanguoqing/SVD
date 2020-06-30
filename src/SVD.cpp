#include "SVD.h"

/*********************************************************************
 * 参数说明：
 * a m*n的实矩阵，返回时其对角线给出奇异值（非递增顺序），其余元素为0
 * m,n 矩阵A的行数和列数
 * u m*m的矩阵，存放左奇异向量
 * v n*n的矩阵，存放右奇异向量
 * eps 双精度实型变量，给定精度要求
 * ka 整形变量，其值为max(m,n)+1
 * 奇异值的情况，此时矩阵A的分解式为UAV，如果返回标志大于0，则说明
 * 程序正常运行
 ********************************************************************/
int SVD(double a[], int m, int n, double u[], double v[], double eps, int ka)
{
	int i, j, k, l, it, ll, kk, ix, iy, mm, nn, iz, ml, ks;
	double d, dd, t, sm, sml, eml, sk, ek, b, c, shh, fg[2], cs[2];
	double *s, *e, *w;
	s = (double *)malloc(ka * sizeof(double));
	e = (double *)malloc(ka * sizeof(double));
	w = (double *)malloc(ka * sizeof(double));
	for (i = 1; i <= m; i++)
	{
		ix = (i - 1) * m + i - 1;
		u[ix] = 0;
	}
	for (i = 1; i <= n; i++)
	{
		iy = (i - 1) * n + i - 1;
		v[iy] = 0;
	}
	it = 60;
	k = n;
	if (m - 1 < n)
		k = m - 1;
	l = m;
	if (n - 2 < m)
		l = n - 2;
	if (l < 0)
		l = 0;
	ll = k;
	if (l > k)
		ll = l;
	/*************** 计算二对角矩阵********************/
	if (ll >= 1)
	{
		for (kk = 1; kk <= ll; kk++) 
		{
			if (kk <= k)
			{
				d = 0.0;
				for (i = kk; i <= m; i++)
				{
					ix = (i - 1) * n + kk - 1; 
					d = d + a[ix] * a[ix]; 
				}
				s[kk - 1] = sqrt(d);  //miu
				if (fabs(s[kk - 1]) > MINNum)
				{
					ix = (kk - 1) * n + kk - 1;
					if (fabs(a[ix]) > MINNum)
					{
						s[kk - 1] = fabs(s[kk - 1]);
						if (a[ix] < 0.0)
							s[kk - 1] = -s[kk - 1];
					}
					for (i = kk; i <= m; i++)
					{
						iy = (i - 1) * n + kk - 1;
						a[iy] = a[iy] / s[kk - 1]; //H1
					}
					a[ix] = 1.0 + a[ix];//这里加1是为了后面计算H1×A方便
				}
				s[kk - 1] = -s[kk - 1];//s为主对角线上的数据
			}
			if (n >= kk + 1)
			{
				for (j = kk + 1; j <= n; j++)
				{
					if ((kk <= k) && (fabs(s[kk - 1]) > MINNum))
					{
						d = 0.0;
						for (i = kk; i <= m; i++)
						{
							ix = (i - 1) * n + kk - 1;
							iy = (i - 1) * n + j - 1;
							d = d + a[ix] * a[iy];
						}
						d = -d / a[(kk - 1) * n + kk - 1];
						for (i = kk; i <= m; i++)
						{
							ix = (i - 1) * n + j - 1;
							iy = (i - 1) * n + kk - 1;
							a[ix] = a[ix] + d * a[iy];
						}
					}
					e[j - 1] = a[(kk - 1) * n + j - 1];// 为副对角线上的数据
				}
			}
			if (kk <= k)
			{
				for (i = kk; i <= m; i++)
				{
					ix = (i - 1) * m + kk - 1;
					iy = (i - 1) * n + kk - 1;
					u[ix] = a[iy];
				}
			}
			if (kk <= l)
			{
				d = 0.0;
				for (i = kk + 1; i <= n; i++)
					d = d + e[i - 1] * e[i - 1];
				e[kk - 1] = sqrt(d);
				if (fabs(e[kk - 1]) > MINNum)
				{
					if (fabs(e[kk]) > MINNum)
					{
						e[kk - 1] = fabs(e[kk - 1]);
						if (e[kk] < 0.0)
							e[kk - 1] = -e[kk - 1];
					}
					for (i = kk + 1; i <= n; i++)
						e[i - 1] = e[i - 1] / e[kk - 1];
					e[kk] = 1.0 + e[kk];
				}
				e[kk - 1] = -e[kk - 1];
				if ((kk + 1 <= m) && (fabs(e[kk - 1]) > MINNum))
				{
					for (i = kk + 1; i <= m; i++)
						w[i - 1] = 0.0;
					for (j = kk + 1; j <= n; j++)
						for (i = kk + 1; i <= m; i++)
							w[i - 1] = w[i - 1] + e[j - 1] * a[(i - 1) * n + j - 1];
					for (j = kk + 1; j <= n; j++)
						for (i = kk + 1; i <= m; i++)
						{
							ix = (i - 1) * n + j - 1;
							a[ix] = a[ix] - w[i - 1] * e[j - 1] / e[kk];
						}
				}
				for (i = kk + 1; i <= n; i++)
					v[(i - 1) * n + kk - 1] = e[i - 1];
			}
		}
	}
	mm = n;
	if (m + 1 < n)
		mm = m + 1;
	if (k < n)
		s[k] = a[k * n + k];
	if (m < mm)
		s[mm - 1] = 0.0;
	if (l + 1 < mm)
		e[l] = a[l * n + mm - 1];
	e[mm - 1] = 0.0;
	nn = m;
	if (m > n)
		nn = n; 
	if (nn >= k + 1)
	{
		for (j = k + 1; j <= nn; j++)
		{
			for (i = 1; i <= m; i++)
				u[(i - 1) * m + j - 1] = 0.0;
			u[(j - 1) * m + j - 1] = 1.0;
		}
	}
	if (k >= 1) 
	{
		for (ll = 1; ll <= k; ll++)
		{
			kk = k - ll + 1; 
			iz = (kk - 1) * m + kk - 1; 
			if (fabs(s[kk - 1]) > MINNum)
			{
				if (nn >= kk + 1)
					for (j = kk + 1; j <= nn; j++)
					{
						d = 0.0;
						for (i = kk; i <= m; i++)
						{
							ix = (i - 1) * m + kk - 1;
							iy = (i - 1) * m + j - 1;
							d = d + u[ix] * u[iy] / u[iz];

						}
						d = -d;
						for (i = kk; i <= m; i++)
						{
							ix = (i - 1) * m + j - 1;
							iy = (i - 1) * m + kk - 1;
							u[ix] = u[ix] + d * u[iy];
						}
					}
				for (i = kk; i <= m; i++)
				{
					ix = (i - 1) * m + kk - 1;
					u[ix] = -u[ix];
				}
				u[iz] = 1.0 + u[iz];
				if (kk - 1 >= 1) 
					for (i = 1; i <= kk - 1; i++)
						u[(i - 1) * m + kk - 1] = 0.0;
			}
			else
			{
				for (i = 1; i <= m; i++)
					u[(i - 1) * m + kk - 1] = 0.0;
				u[(kk - 1) * m + kk - 1] = 1.0;
			}
		}
	}
	/*************** 计算二对角矩阵********************/
	for (ll = 1; ll <= n; ll++)
	{
		kk = n - ll + 1;
		iz = kk * n + kk - 1;
		//if((kk<=l)&&(e[kk-1]!=0.0))/////////////////////////////
		if ((kk <= l) && (fabs(e[kk - 1]) > MINNum))
		{
			for (j = kk + 1; j <= n; j++)
			{
				d = 0.0;
				for (i = kk + 1; i <= n; i++)
				{
					ix = (i - 1) * n + kk - 1;
					iy = (i - 1) * n + j - 1;
					d = d + v[ix] * v[iy] / v[iz];
				}
				d = -d;
				for (i = kk + 1; i <= n; i++)
				{
					ix = (i - 1) * n + j - 1;
					iy = (i - 1) * n + kk - 1;
					v[ix] = v[ix] + d * v[iy];
				}
			}
		}
		for (i = 1; i <= n; i++)
			v[(i - 1) * n + kk - 1] = 0.0;
		v[iz - n] = 1.0;
	}
	for (i = 1; i <= m; i++)
		for (j = 1; j <= n; j++)
			a[(i - 1) * n + j - 1] = 0.0;
	ml = mm;
	it = 60;
	while (1 == 1) //////////////////////////////////
	{
		if (mm == 0)
		{
			QRCompute(a, e, s, v, m, n);
			free(s);
			free(e);
			free(w);
			return l;
		}
		if (it == 0)
		{
			QRCompute(a, e, s, v, m, n);
			free(s);
			free(e);
			free(w);
			return -1;
		}
		kk = mm - 1;
		// 收敛判断,置0
		while ((kk != 0) && (fabs(e[kk - 1]) > MINNum))
		{
			d = fabs(s[kk - 1]) + fabs(s[kk]);
			dd = fabs(e[kk - 1]);
			if (dd > eps * d)
				kk = kk - 1;
			else
				e[kk - 1] = 0.0;
		}
		if (kk == mm - 1)
		{
			kk = kk + 1;
			if (s[kk - 1] < 0.0)
			{
				s[kk - 1] = -s[kk - 1];
				for (i = 1; i <= n; i++)
				{
					ix = (i - 1) * n + kk - 1;
					v[ix] = -v[ix];
				}
			}
			while ((kk != ml) && (s[kk - 1] < s[kk]))
			{
				d = s[kk - 1];
				s[kk - 1] = s[kk];
				s[kk] = d;
				if (kk < n)
					for (i = 1; i <= n; i++)
					{
						ix = (i - 1) * n + kk - 1;
						iy = (i - 1) * n + kk;
						d = v[ix];
						v[ix] = v[iy];
						v[iy] = d;
					}
				if (kk < m)
					for (i = 1; i <= m; i++)
					{
						ix = (i - 1) * m + kk - 1;
						iy = (i - 1) * m + kk;
						d = u[ix];
						u[ix] = u[iy];
						u[iy] = d;
					}
				kk = kk + 1;
			}
			it = 60;
			mm = mm - 1;
		}
		else
		{
			ks = mm;
			while ((ks > kk) && (fabs(s[ks - 1]) > MINNum))
			{
				d = 0.0;
				if (ks != mm)
					d = d + fabs(e[ks - 1]);
				if (ks != kk + 1)
					d = d + fabs(e[ks - 2]);
				dd = fabs(s[ks - 1]);
				if (dd > eps * d)
					ks = ks - 1;
				else
					s[ks - 1] = 0.0;
			}
			if (ks == kk)
			{
				kk = kk + 1;
				d = fabs(s[mm - 1]);
				t = fabs(s[mm - 2]);
				if (t > d)
					d = t;
				t = fabs(e[mm - 2]);
				if (t > d)
					d = t;
				t = fabs(s[kk - 1]);
				if (t > d)
					d = t;
				t = fabs(e[kk - 1]);
				if (t > d)
					d = t;
				// 计算QR分解的位移
				sm = s[mm - 1] / d;
				sml = s[mm - 2] / d;
				eml = e[mm - 2] / d;
				sk = s[kk - 1] / d;
				ek = e[kk - 1] / d;
				b = ((sml + sm) * (sml - sm) + eml * eml) / 2.0;
				c = sm * eml;
				c = c * c;
				shh = 0.0;
				if ((fabs(b) > MINNum) || (fabs(c) > MINNum))
				{
					shh = sqrt(b * b + c);
					if (b < 0.0)
						shh = -shh;
					shh = c / (b + shh);
				}
				fg[0] = (sk + sm) * (sk - sm) - shh;//x
				fg[1] = sk * ek;//y
				for (i = kk; i <= mm - 1; i++)
				{
					Givens(fg, cs); //计算Givens变换矩阵参数
					if (i != kk)
						e[i - 2] = fg[0];
					fg[0] = cs[0] * s[i - 1] + cs[1] * e[i - 1];
					e[i - 1] = cs[0] * e[i - 1] - cs[1] * s[i - 1];
					fg[1] = cs[1] * s[i];
					s[i] = cs[0] * s[i];
					if ((fabs(cs[0] - 1.0) > MINNum) || (fabs(cs[1]) > MINNum)) // s c最好不为零
						for (j = 1; j <= n; j++)
						{
							ix = (j - 1) * n + i - 1;
							iy = (j - 1) * n + i;
							d = cs[0] * v[ix] + cs[1] * v[iy];
							v[iy] = -cs[1] * v[ix] + cs[0] * v[iy];
							v[ix] = d;
						}
					Givens(fg, cs);
					s[i - 1] = fg[0];
					fg[0] = cs[0] * e[i - 1] + cs[1] * s[i];
					s[i] = -cs[1] * e[i - 1] + cs[0] * s[i];
					fg[1] = cs[1] * e[i];
					e[i] = cs[0] * e[i];
					if (i < m)
						if ((fabs(cs[0] - 1.0) > MINNum) || (fabs(cs[1]) > MINNum))
							for (j = 1; j <= m; j++)
							{
								ix = (j - 1) * m + i - 1;
								iy = (j - 1) * m + i;
								d = cs[0] * u[ix] + cs[1] * u[iy];
								u[iy] = -cs[1] * u[ix] + cs[0] * u[iy];
								u[ix] = d;
							}
				}
				e[mm - 2] = fg[0];
				it = it - 1;
			}
			else
			{
				if (ks == mm)
				{
					kk = kk + 1;
					fg[1] = e[mm - 2];
					e[mm - 2] = 0.0;
					for (ll = kk; ll <= mm - 1; ll++)
					{
						i = mm + kk - ll - 1;
						fg[0] = s[i - 1];
						Givens(fg, cs);
						s[i - 1] = fg[0];
						if (i != kk)
						{
							fg[1] = -cs[1] * e[i - 2];
							e[i - 2] = cs[0] * e[i - 2];
						}
						if ((fabs(cs[0] - 1.0) > MINNum) || (fabs(cs[1]) > MINNum))
							for (j = 1; j <= n; j++)
							{
								ix = (j - 1) * n + i - 1;
								iy = (j - 1) * n + mm - 1;
								d = cs[0] * v[ix] + cs[1] * v[iy];
								v[iy] = -cs[1] * v[ix] + cs[0] * v[iy];
								v[ix] = d;
							}
					}
				}
				else
				{
					kk = ks + 1;
					fg[1] = e[kk - 2];
					e[kk - 2] = 0.0;
					for (i = kk; i <= mm; i++)
					{
						fg[0] = s[i - 1];
						Givens(fg, cs);
						s[i - 1] = fg[0];
						fg[1] = -cs[1] * e[i - 1];
						e[i - 1] = cs[0] * e[i - 1];
						if ((fabs(cs[0] - 1.0) > MINNum) || (fabs(cs[1]) > MINNum))
							for (j = 1; j <= m; j++)
							{
								ix = (j - 1) * m + i - 1;
								iy = (j - 1) * m + kk - 2;
								d = cs[0] * u[ix] + cs[1] * u[iy];
								u[iy] = -cs[1] * u[ix] + cs[0] * u[iy];
								u[ix] = d;
							}
					}
				}
			}
		}
	}
	free(s);
	free(e);
	free(w);
	return l;
}
void QRCompute(double a[], double e[], double s[], double v[], int m, int n)
{
	int i, j, p, q;
	double d;
	if (m >= n)
		i = n;
	else
		i = m;
	for (j = 1; j <= i - 1; j++)
	{
		a[(j - 1) * n + j - 1] = s[j - 1];
		a[(j - 1) * n + j] = e[j - 1];
	}
	a[(i - 1) * n + i - 1] = s[i - 1];
	if (m < n)
		a[(i - 1) * n + i] = e[i - 1];
	for (i = 1; i <= n - 1; i++)
		for (j = i + 1; j <= n; j++)
		{
			p = (i - 1) * n + j - 1;
			q = (j - 1) * n + i - 1;
			d = v[p];
			v[p] = v[q];
			v[q] = d;
		}
	return;
}
void Givens(double fg[2], double cs[2])
{
	double r, d;
	if ((fabs(fg[0]) + fabs(fg[1])) < MINNum)
	{
		cs[0] = 1.0;
		cs[1] = 0.0;
		d = 0.0;
	}
	else
	{
		d = sqrt(fg[0] * fg[0] + fg[1] * fg[1]);
		if (fabs(fg[0]) > fabs(fg[1]))
		{
			d = fabs(d);
			if (fg[0] < 0.0)
				d = -d;
		}
		if (fabs(fg[1]) >= fabs(fg[0]))
		{
			d = fabs(d);
			if (fg[1] < 0.0)
				d = -d;
		}
		cs[0] = fg[0] / d;
		cs[1] = fg[1] / d;
	}
	r = 1.0;
	if (fabs(fg[0]) > fabs(fg[1]))
		r = cs[1];
	else
		if (fabs(cs[0]) > MINNum)
		r = 1.0 / cs[0];
	fg[0] = d;
	fg[1] = r;
	return;
}
void TriMul(double a[], double b[], int m, int n, int k, double c[])
{
	int i, j, l, u;
	for (i = 0; i <= m - 1; i++)
		for (j = 0; j <= k - 1; j++)
		{
			u = i * k + j;
			c[u] = 0;
			for (l = 0; l <= n - 1; l++)
				c[u] = c[u] + a[i * n + l] * b[l * k + j];
		}
	return;
}
