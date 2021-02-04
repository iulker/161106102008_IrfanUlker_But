#include <iostream>
#include <math.h>
#include <conio.h>
#include <fstream>
using namespace std;


void matrisYazdir(double** A, int i1, int i2)
{
	for (int i = 0; i < i1; ++i)
	{
		for (int j = 0; j < i2; ++j)
		{
			cout << A[i][j] << "\t";
		}
		cout << endl;
	}
}

void matrisYazdir(double* bT, int n)
{
	for (int i = 0; i < n; i++)
	{
		cout << bT[i] << "\t";
	}
}
// 1- A matrisinin determinantýný hesaplayýnýz.
double determinantHesapla(double** A, int n)
{
	double det = 0;
	double** submatrix = new double* [n];
	for (int i = 0; i < n; i++) {
		submatrix[i] = new double[n];
	}
	if (n == 2)
		return ((A[0][0] * A[1][1]) - (A[1][0] * A[0][1]));
	else {
		for (int x = 0; x < n; x++) {
			int subi = 0;
			for (int i = 1; i < n; i++) {
				int subj = 0;
				for (int j = 0; j < n; j++) {
					if (j == x)
						continue;
					submatrix[subi][subj] = A[i][j];
					subj++;
				}
				subi++;
			}
			det = det + (pow(-1, x) * A[0][x] * determinantHesapla(submatrix, n - 1));
		}
	}
	return det;
}

//2- ATA matrisinin izini hesaplayýnýz.
double matrisinIziniHesaplama(double** A, int i1, int i2)
{
	double** devrik = new double* [i1];
	for (int i = 0; i < i1; i++) {
		devrik[i] = new double[i2];
	}
	double** carpim = new double* [i1];
	for (int i = 0; i < i1; i++) {
		carpim[i] = new double[i2];
	}
	double toplam = 0;
	double ara = 0;
	for (int i = 0; i < i1; i++)
	{
		for (int j = 0; j < i2; j++)
		{
			devrik[i][j] = A[j][i];
		}
	}
	//ATA matrisini buluyoruz
	for (int i = 0; i < i1; i++)
	{
		for (int j = 0; j < i2; j++)
		{
			for (int k = 0; k < i2; k++)
			{
				ara += (devrik[i][k] * A[k][j]);
			}
			carpim[i][j] = ara;
			ara = 0;
		}
	}
	//ATA matrisinin izini buluyoruz
	for (int i = 0; i < i1; i++)
	{
		for (int j = 0; j < i2; j++)
		{
			if (i == j)
			{
				toplam += carpim[i][j];
			}
		}
	}
	return toplam;
}
//3- A matrisinin satýr normlarýný hesaplayýnýz.
void satirNormlari(double** A, int i1, int i2)
{
	double norm = 0;
	for (int i = 0; i < i1; i++)
	{
		for (int j = 0; j < i2; j++)
		{
			norm += A[i][j];
		}
		cout << i + 1 << ". satirin normu -->" << norm << endl;
		norm = 0;
	}
}
//4- A matrisinin sutun normlarýný hesaplayýnýz.
void sutunNormlari(double** A, int i1, int i2)
{
	double norm = 0;
	for (int i = 0; i < i2; i++)
	{
		for (int j = 0; j < i1; j++)
		{
			norm += A[j][i];
		}
		cout << i + 1 << ". sutunun normu -->" << norm << endl;
		norm = 0;
	}
}
//5-) A matrisinin Öklid normunu (N(A)) hesaplayýnýz.
double oklidNormu(double** A, int i1, int i2)
{
	double norm = 0, aradeger = 0;
	for (int i = 0; i < i1; i++)
	{
		for (int j = 0; j < i2; j++)
		{
			aradeger = pow(A[i][j], 2);
			norm += aradeger;
		}
	}
	return sqrt(norm);
}
//6-) N(A)=(iz(ATA))^(1/2) olduðunun saðlamasýný gerçekleþtiriniz.
void saglama(double** A, int i1, int i2)
{
	double oklid = oklidNormu(A, i1, i2);
	double iz = matrisinIziniHesaplama(A, i1, i2);
	cout << "6-) N(A)=(iz(ATA))^(1/2) cunku;" << endl;
	cout << " " << oklid << "=" << iz << "^(1/2)" << endl;
	cout << "Yukaridaki islemi yaptigimizda --> " << oklid << "=" << sqrt(iz) << endl;
}
//7-) A matrisini Öklid normuna göre normlaþtýrýnýz.
void normlastirma(double** A, int i1, int i2)
{
	double** Anormal = new double* [i1];

	for (int i = 0; i < i1; i++) {
		Anormal[i] = new double[i2];
	}

	for (int i = 0; i < i1; i++)
	{
		for (int j = 0; j < i2; j++)
		{
			Anormal[i][j] = A[i][j] / oklidNormu(A, i1, i2);
		}
	}
	matrisYazdir(Anormal, i1, i2);
	cout << endl;
}
//8-) A matrisinin özdeðerlerini hesaplayýnýz.
//9-) A matrisinin Spektral (Todd) þart sayýsýný hesaplayarak kararsýzlýðýný yorumlayýnýz.
//10-) A matrisinin Hadamard þart sayýsýný hesaplayarak kararsýzlýðýný yorumlayýnýz
void hadamard(double** A, int n)
{
	double d = determinantHesapla(A, n);
	double* alpha = new double[n];
	double aratoplam = 0, carpim = 1;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			aratoplam += pow(A[i][j], 2);
		}
		alpha[i] = sqrt(aratoplam);
		aratoplam = 0;
	}
	for (int i = 0; i < n; i++)
	{
		carpim *= alpha[i];
	}
	double hadamard = d / carpim;
	cout << "10-) A matrisinin Hadamard sart sayisi -->" << hadamard;
	if (hadamard < 0.01)
	{
		cout << " ve kararsizdir " << endl;
	}
	else
	{
		cout << " ve kararlidir" << endl;
	}
}
//11-) Kramer kuralý ile A matrisinin tersini hesaplayýnýz.
void kramerTers(double** A, int n)
{
	double anadet = determinantHesapla(A, n);
	double** tersmatris = new double* [n];
	for (int i = 0; i < n; i++) {
		tersmatris[i] = new double[n];
	}
	for (int x = 0; x < n; x++) {
		for (int y = 0; y < n; y++) {


			int n1 = n - 1;
			double** araMatris = new double* [n1];
			for (int a = 0; a < n1; a++) {
				araMatris[a] = new double[n1];
			}

			int i1 = 0;
			int i2 = 0;

			for (int i = 0; i < n; i++) {
				if (i == x) continue;
				for (int j = 0; j < n; j++) {
					if (j == y) continue;

					araMatris[i1][i2] = A[i][j];

					i2++;
				}
				i1++;
			}

			double aradet = determinantHesapla((double**)araMatris, n - 1);

			tersmatris[x][y] = aradet / anadet;

		}
	}

	//double aradet = 0;
	//double araMatris1[4][4] = { {A[1][1] ,A[1][2] ,A[1][3] ,A[1][4]},
	//{ A[2][1] ,A[2][2] ,A[2][3] ,A[2][4]},
	//{ A[3][1] ,A[3][2] ,A[3][3] ,A[3][4]},
	//{ A[4][1] ,A[4][2] ,A[4][3] ,A[4][4]} };
	//aradet = determinantHesapla((double**)araMatris1, n);
	//double c11 = aradet;
	//double araMatris2[4][4] = { {A[1][0] ,A[1][2] ,A[1][3] ,A[1][4]},
	//{ A[2][0] ,A[2][2] ,A[2][3] ,A[2][4]},
	//{ A[3][0] ,A[3][2] ,A[3][3] ,A[3][4]},
	//{ A[4][0] ,A[4][2] ,A[4][3] ,A[4][4]} };
	//aradet = determinantHesapla((double**)araMatris2, n);
	//double c12 = -aradet;
	//double araMatris3[4][4] = { {A[1][0] ,A[1][1] ,A[1][3] ,A[1][4]},
	//{ A[2][0] ,A[2][1] ,A[2][3] ,A[2][4]},
	//{ A[3][0] ,A[3][1] ,A[3][3] ,A[3][4]},
	//{ A[4][0] ,A[4][1] ,A[4][3] ,A[4][4]} };
	//aradet = determinantHesapla4x4(araMatris3, 4);
	//double c13 = aradet;
	//double araMatris4[4][4] = { {A[1][0] ,A[1][1] ,A[1][2] ,A[1][4]},
	//{ A[2][0] ,A[2][1] ,A[2][2] ,A[2][4]},
	//{ A[3][0] ,A[3][1] ,A[3][2] ,A[3][4]},
	//{ A[4][0] ,A[4][1] ,A[4][2] ,A[4][4]} };
	//aradet = determinantHesapla4x4(araMatris4, 4);
	//double c14 = -aradet;
	//double araMatris5[4][4] = { {A[1][0] ,A[1][1] ,A[1][2] ,A[1][3]},
	//{ A[2][0] ,A[2][1] ,A[2][2] ,A[2][3]},
	//{ A[3][0] ,A[3][1] ,A[3][2] ,A[3][3]},
	//{ A[4][0] ,A[4][1] ,A[4][2] ,A[4][3]} };
	//aradet = determinantHesapla4x4(araMatris5, 4);
	//double c15 = aradet;
	//double araMatris6[4][4] = { {A[0][1] ,A[0][2] ,A[0][3] ,A[0][4]},
	//{ A[2][1] ,A[2][2] ,A[2][3] ,A[2][4]},
	//{ A[3][1] ,A[3][2] ,A[3][3] ,A[3][4]},
	//{ A[4][1] ,A[4][2] ,A[4][3] ,A[4][4]} };
	//aradet = determinantHesapla4x4(araMatris6, 4);
	//double c21 = -aradet;
	//double araMatris7[4][4] = { {A[0][0] ,A[0][2] ,A[0][3] ,A[0][4]},
	//{ A[2][0] ,A[2][2] ,A[2][3] ,A[2][4]},
	//{ A[3][0] ,A[3][2] ,A[3][3] ,A[3][4]},
	//{ A[4][0] ,A[4][2] ,A[4][3] ,A[4][4]} };
	//aradet = determinantHesapla4x4(araMatris7, 4);
	//double c22 = aradet;
	//double araMatris8[4][4] = { {A[0][0] ,A[0][1] ,A[0][3] ,A[0][4]},
	//{ A[2][0] ,A[2][1] ,A[2][3] ,A[2][4]},
	//{ A[3][0] ,A[3][1] ,A[3][3] ,A[3][4]},
	//{ A[4][0] ,A[4][1] ,A[4][3] ,A[4][4]} };
	//aradet = determinantHesapla4x4(araMatris8, 4);
	//double c23 = -aradet;
	//double araMatris9[4][4] = { {A[0][0] ,A[0][1] ,A[0][2] ,A[0][4]},
	//{ A[2][0] ,A[2][1] ,A[2][2] ,A[2][4]},
	//{ A[3][0] ,A[3][1] ,A[3][2] ,A[3][4]},
	//{ A[4][0] ,A[4][1] ,A[4][2] ,A[4][4]} };
	//aradet = determinantHesapla4x4(araMatris9, 4);
	//double c24 = aradet;
	//double araMatris10[4][4] = { {A[0][0] ,A[0][1] ,A[0][2] ,A[0][3]},
	//{ A[2][0] ,A[2][1] ,A[2][2] ,A[2][3]},
	//{ A[3][0] ,A[3][1] ,A[3][2] ,A[3][3]},
	//{ A[4][0] ,A[4][1] ,A[4][2] ,A[4][3]} };
	//aradet = determinantHesapla4x4(araMatris10, 4);
	//double c25 = -aradet;
	//double araMatris11[4][4] = { {A[0][1] ,A[0][2] ,A[0][3] ,A[0][4]},
	//{ A[1][1] ,A[1][2] ,A[1][3] ,A[1][4]},
	//{ A[3][1] ,A[3][2] ,A[3][3] ,A[3][4]},
	//{ A[4][1] ,A[4][2] ,A[4][3] ,A[4][4]} };
	//aradet = determinantHesapla4x4(araMatris11, 4);
	//double c31 = aradet;
	//double araMatris12[4][4] = { {A[0][0] ,A[0][2] ,A[0][3] ,A[0][4]},
	//{ A[1][0] ,A[1][2] ,A[1][3] ,A[1][4]},
	//{ A[3][0] ,A[3][2] ,A[3][3] ,A[3][4]},
	//{ A[4][0] ,A[4][2] ,A[4][3] ,A[4][4]} };
	//aradet = determinantHesapla4x4(araMatris12, 4);
	//double c32 = -aradet;
	//double araMatris13[4][4] = { {A[0][0] ,A[0][1] ,A[0][3] ,A[0][4]},
	//{ A[1][0] ,A[1][1] ,A[1][3] ,A[1][4]},
	//{ A[3][0] ,A[3][1] ,A[3][3] ,A[3][4]},
	//{ A[4][0] ,A[4][1] ,A[4][3] ,A[4][4]} };
	//aradet = determinantHesapla4x4(araMatris13, 4);
	//double c33 = aradet;
	//double araMatris14[4][4] = { {A[0][0] ,A[0][1] ,A[0][2] ,A[0][4]},
	//{ A[1][0] ,A[1][1] ,A[1][2] ,A[1][4]},
	//{ A[3][0] ,A[3][1] ,A[3][2] ,A[3][4]},
	//{ A[4][0] ,A[4][1] ,A[4][2] ,A[4][4]} };
	//aradet = determinantHesapla4x4(araMatris14, 4);
	//double c34 = -aradet;
	//double araMatris15[4][4] = { {A[0][0] ,A[0][1] ,A[0][2] ,A[0][3]},
	//{ A[1][0] ,A[1][1] ,A[1][2] ,A[1][3]},
	//{ A[3][0] ,A[3][1] ,A[3][2] ,A[3][3]},
	//{ A[4][0] ,A[4][1] ,A[4][2] ,A[4][3]} };
	//aradet = determinantHesapla4x4(araMatris15, 4);
	//double c35 = aradet;
	//double araMatris16[4][4] = { {A[0][1] ,A[0][2] ,A[0][3] ,A[0][4]},
	//{ A[1][1] ,A[1][2] ,A[1][3] ,A[1][4]},
	//{ A[2][1] ,A[2][2] ,A[2][3] ,A[2][4]},
	//{ A[4][1] ,A[4][2] ,A[4][3] ,A[4][4]} };
	//aradet = determinantHesapla4x4(araMatris16, 4);
	//double c41 = -aradet;
	//double araMatris17[4][4] = { {A[0][0] ,A[0][2] ,A[0][3] ,A[0][4]},
	//{ A[1][0] ,A[1][2] ,A[1][3] ,A[1][4]},
	//{ A[2][0] ,A[2][2] ,A[2][3] ,A[2][4]},
	//{ A[4][0] ,A[4][2] ,A[4][3] ,A[4][4]} };
	//aradet = determinantHesapla4x4(araMatris17, 4);
	//double c42 = aradet;
	//double araMatris18[4][4] = { {A[0][0] ,A[0][1] ,A[0][3] ,A[0][4]},
	//{ A[1][0] ,A[1][1] ,A[1][3] ,A[1][4]},
	//{ A[2][0] ,A[2][1] ,A[2][3] ,A[2][4]},
	//{ A[4][0] ,A[4][1] ,A[4][3] ,A[4][4]} };
	//aradet = determinantHesapla4x4(araMatris18, 4);
	//double c43 = -aradet;
	//double araMatris19[4][4] = { {A[0][0] ,A[0][1] ,A[0][2] ,A[0][4]},
	//{ A[1][0] ,A[1][1] ,A[1][2] ,A[1][4]},
	//{ A[2][0] ,A[2][1] ,A[2][2] ,A[2][4]},
	//{ A[4][0] ,A[4][1] ,A[4][2] ,A[4][4]} };
	//aradet = determinantHesapla4x4(araMatris19, 4);
	//double c44 = aradet;
	//double araMatris20[4][4] = { {A[0][0] ,A[0][1] ,A[0][2] ,A[0][3]},
	//{ A[1][0] ,A[1][1] ,A[1][2] ,A[1][3]},
	//{ A[2][0] ,A[2][1] ,A[2][2] ,A[2][3]},
	//{ A[4][0] ,A[4][1] ,A[4][2] ,A[4][3]} };
	//aradet = determinantHesapla4x4(araMatris20, 4);
	//double c45 = -aradet;
	//double araMatris21[4][4] = { {A[0][1] ,A[0][2] ,A[0][3] ,A[0][4]},
	//{ A[1][1] ,A[1][2] ,A[1][3] ,A[1][4]},
	//{ A[2][1] ,A[2][2] ,A[2][3] ,A[2][4]},
	//{ A[3][1] ,A[3][2] ,A[3][3] ,A[3][4]} };
	//aradet = determinantHesapla4x4(araMatris21, 4);
	//double c51 = aradet;
	//double araMatris22[4][4] = { {A[0][0] ,A[0][2] ,A[0][3] ,A[0][4]},
	//{ A[1][0] ,A[1][2] ,A[1][3] ,A[1][4]},
	//{ A[2][0] ,A[2][2] ,A[2][3] ,A[2][4]},
	//{ A[3][0] ,A[3][2] ,A[3][3] ,A[3][4]} };
	//aradet = determinantHesapla4x4(araMatris22, 4);
	//double c52 = -aradet;
	//double araMatris23[4][4] = { {A[0][0] ,A[0][1] ,A[0][3] ,A[0][4]},
	//{ A[1][0] ,A[1][1] ,A[1][3] ,A[1][4]},
	//{ A[2][0] ,A[2][1] ,A[2][3] ,A[2][4]},
	//{ A[3][0] ,A[3][1] ,A[3][3] ,A[3][4]} };
	//aradet = determinantHesapla4x4(araMatris23, 4);
	//double c53 = aradet;
	//double araMatris24[4][4] = { {A[0][0] ,A[0][1] ,A[0][2] ,A[0][4]},
	//{ A[1][0] ,A[1][1] ,A[1][2] ,A[1][4]},
	//{ A[2][0] ,A[2][1] ,A[2][2] ,A[2][4]},
	//{ A[3][0] ,A[3][1] ,A[3][2] ,A[3][4]} };
	//aradet = determinantHesapla4x4(araMatris24, 4);
	//double c54 = -aradet;
	//double araMatris25[4][4] = { {A[0][0] ,A[0][1] ,A[0][2] ,A[0][3]},
	//{ A[1][0] ,A[1][1] ,A[1][2] ,A[1][3]},
	//{ A[2][0] ,A[2][1] ,A[2][2] ,A[2][3]},
	//{ A[3][0] ,A[3][1] ,A[3][2] ,A[3][3]} };
	//aradet = determinantHesapla4x4(araMatris25, 4);
	//double c55 = aradet;

	/*double tersmatris[5][5] = { {c11 / anadet ,c21 / anadet ,c31 / anadet ,c41 / anadet ,c51 / anadet},
	{c12 / anadet ,c22 / anadet ,c32 / anadet ,c42 / anadet ,c52 / anadet},
	{c13 / anadet ,c23 / anadet ,c33 / anadet ,c43 / anadet ,c53 / anadet},
	{c14 / anadet ,c24 / anadet ,c34 / anadet ,c44 / anadet ,c54 / anadet},
	{c15 / anadet ,c25 / anadet ,c35 / anadet ,c45 / anadet ,c55 / anadet} };*/
	matrisYazdir(tersmatris, n, n);
}
//14-) Gauss algoritmasý ile x bilinmeyenler vektörünü hesaplayýnýz. Önce alt üçgen matris oluþturulur.
void gaussIleXBilinmeyenlerVektoru(double** B, int i1, int i2)
{
	double sifir = 0.0;
	double carpim = 0.0;
	for (int k = 0; k <= i1 - 1; k++)
	{
		for (int i = k + 1; i <= i1 - 1; i++)
		{
			sifir = B[i][k] / B[k][k];
			for (int j = k; j <= i2 - 1; j++)
			{
				carpim = B[k][j] * sifir;
				B[i][j] = B[i][j] - carpim;
			}
			carpim = 0.0;
		}
		sifir = 0.0;
	}
	B[3][1] = 0.0;
	double* x = new double[i2 - 1];
	x[i2 - 2] = B[i1 - 1][i2 - 1] / B[i1 - 1][i2 - 2];
	for (int i = i1 - 2; i >= 0; i--)
	{
		double toplam = 0;
		for (int j = i + 1; j <= i2 - 2; j++)
		{
			toplam = B[i][j] * x[j];
		}
		double tmp = B[i][i2 - 1] - toplam;
		x[i] = tmp / B[i][i];
	}
	matrisYazdir(x, i2 - 1);
}
//15-) Gauss Jordan
void gaussJordan(double** B, int i1, int i2)
{
	double b;//Line 1
	int i, j, k;
	double* x = new double[i1];
	for (j = 0; j <= i1 - 1; j++)
	{
		for (i = 0; i <= i1 - 1; i++)
		{
			if (i != j)
			{
				b = B[i][j] / B[j][j];
				for (k = 0; k <= i1; k++)
				{
					B[i][k] = B[i][k] - b * B[j][k];
				}
			}
		}
	}
	for (i = 0; i <= i1 - 1; i++)
	{
		x[i] = B[i][i2 - 1] / B[i][i];
		cout << x[i] << "\t";
	}
}
int main()
{
	bool isExit = false;
	int secim = 0;
	int m;
	double** A;
	//double** B = new double* [m];

	while (secim != 1 && secim != 2) {
		cout << "Matris girisi nasil olsun:" << endl;
		cout << "	1) Dosyadan oku(A.txt):" << endl;
		cout << "	2) Elle gir:" << endl;
		cin >> secim;
	}

	if (secim == 1) {
		ifstream f("A.txt");
		f >> m;
		A = new double* [m];
		for (int i = 0; i < m; i++) {
			A[i] = new double[m];
		}
		for (int i = 0; i < m; i++)
			for (int j = 0; j < m; j++)
				f >> A[i][j];
		cout << "A.txt dosyasindan okundu"<<endl;
		/*ifstream f2("B.txt");
		f2 >> m >> n;
		for (int i = 0; i < m; i++) {
			B[i] = new double[n];
		}
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
				f2 >> B[i][j];*/
	}
	else if (secim == 2)
	{
		cout << "Kare matris icin boyut girin:" << endl;
		cin >> m;
		A = new double* [m];
		for (int i = 0; i < m; i++) {
			A[i] = new double[m];
		}
		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < m; j++)
			{
				cout << "A[" << i << "][" << j << "] icin deger girin:";
				cin >> A[i][j];
			}
		}
	}
	else
	{
		return -1;
	}
	do {

		system("cls");

		cout << "A matrisi:" << endl;

		matrisYazdir(A, m, m);

		//cout << endl << "B.txt dosyasindan okunan B matrisi:" << endl;
		//matrisYazdir(B, m, n);
		//cout << endl;
		/*double** A = new double* [] {
			new double[] {0.5470, 0.1835, 0.9294, 0.3063, 0.6443},
				new double[] {0.2963, 0.3685, 0.7757, 0.5085, 0.3786},
				new double[] {0.7447, 0.6256, 0.4868, 0.5108, 0.8116},
				new double[] {0.1890, 0.7802, 0.4359, 0.8176, 0.5328},
				new double[] {0.6868, 0.0811, 0.4468, 0.7948, 0.3507}
		};

		double* bT = new double[] { 0.9390, 0.8759, 0.5502, 0.6225, 0.5870 };

		double** B = new double* [] {
			new double[] {0.5470, 0.1835, 0.9294, 0.3063, 0.6443, 0.9390},
				new double[] {0.2963, 0.3685, 0.7757, 0.5085, 0.3786, 0.8759},
				new double[] {0.7447, 0.6256, 0.4868, 0.5108, 0.8116, 0.5502},
				new double[] {0.1890, 0.7802, 0.4359, 0.8176, 0.5328, 0.6225},
				new double[] {0.6868, 0.0811, 0.4468, 0.7948, 0.3507, 0.5870}
		};*/

		cout << "Menu:" << endl;
		cout << "1 - A matrisinin determinanti" << endl;
		cout << "2 - A matrisinin transpozunun izi" << endl;
		cout << "3 - A matrisinin satir normlari" << endl;
		cout << "4 - A matrisinin sutun normlari" << endl;
		cout << "5 - A matrisinin oklid normu" << endl;
		cout << "6 - A matrisinin determinanti" << endl;
		cout << "7 - A matrisinin oklid normu ile normallastirilmis hali" << endl;
		//cout << "8 - A matrisinin determinanti" << endl;
		//cout << "9 - A matrisinin determinanti" << endl;
		cout << "10 - A matrisinin Hadamard sart sayisi ve kararliligi" << endl;
		cout << "11 - A matrisini kramer kurali ile tersi hesapla" << endl;
		//cout << "12 - A matrisinin determinantý" << endl;
		//cout << "13 - A matrisinin determinantý" << endl;
		cout << "14 - Gauss algoritmasi ile x bilinmeyenler vektoru hesapla" << endl;
		cout << "15 - Gauss jordan yontemi ile x bilinmeyenler vektoru hesapla" << endl;
		cout << "0 - Exit" << endl;
		cout << "Lutfen seciminizi girin(1-15): ";

		cin >> secim;
		system("cls");

		switch (secim)
		{
		case 1:
			cout << "1-) A matrisinin determinanti --> " << determinantHesapla((double**)A, m) << endl;
			break;
		case 2:
			cout << "2-) A matrisinin transpozunun izi -->" << matrisinIziniHesaplama((double**)A, m, m) << endl;
			break;
		case 3:
			cout << "3-) A matrisinin satir normlari ;" << endl;
			satirNormlari((double**)A, m, m);
			break;
		case 4:
			cout << "4-) A matrisinin sutun normlari ;" << endl;
			sutunNormlari((double**)A, m, m);
			break;
		case 5:
			cout << "5-) A matrisinin oklid normu --> " << oklidNormu((double**)A, m, m) << endl;
		case 6:
			saglama((double**)A, m, m);
			break;
		case 7:
			cout << "7-) A matrisinin oklid normu ile normallastirilmis hali;" << endl;
			normlastirma((double**)A, m, m);
			break;
		case 8:
			break;
		case 9:
			break;
		case 10:
			hadamard((double**)A, m);
			break;
		case 11:
			cout << "11-) A matrisinin kramer kurali ile tersi hesaplandi ve matris;" << endl;
			kramerTers((double**)A, m);
			break;
		case 12:
			break;
		case 13:
			break;
		case 14:
			cout << "14-)Gauss algoritmasi ile x bilinmeyenler vektoru hesaplandi ve matris;" << endl;
			gaussIleXBilinmeyenlerVektoru((double**)A, m, m);
			break;
		case 15:
			cout << "15-)Gauss jordan yontemi ile x bilinmeyenler vektoru hesaplandi ve matris;" << endl;
			gaussJordan((double**)A, m, m);
			break;
		case 0:
			isExit = true;
			break;
		default:
			cout << "Hatali Secim Yaptiniz" << endl;
			break;
		}

		cout << endl; // << endl << "Devam etmek için bir tuþa basýn" << endl;

		if (!isExit) system("pause");

	} while (!isExit);

	return 0;
}