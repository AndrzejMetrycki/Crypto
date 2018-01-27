#include <stdint.h>
#include <iostream>

#include <thread>
#include <vector>
#include <mutex>
#include <chrono>
#include <cstdint>
#include <math.h>

#include <valarray>


using namespace std;
using namespace chrono;

//__int128 mint; 
typedef uint64_t _int_t;

class CPrimeNr_Virtual
{
public:
	virtual bool IsPrime(_int_t iNr) = 0;
};




class CPrimeNr : public CPrimeNr_Virtual
{
public:


	virtual void OnPrimeNumberFound(_int_t iNr)
	{
		static mutex mMutex;
		static int iCounter = 0;

		mMutex.lock();
		cout << "Liczba Pierwsza[" << iCounter++ << "]: " << iNr << endl;
		mMutex.unlock();

	}

	void FindPrimeNumbers(_int_t iBegin, _int_t iEnd)
	{
		for (_int_t i = iBegin; i <= iEnd; ++i)
		{
			if (IsPrime(i))
			{
				OnPrimeNumberFound(i);
			}
		}
	}



	void Multithread_FindPrimaryNumbers(_int_t iBegin, _int_t iEnd)
	{
		vector<thread> vThd;
		uint64_t iTotalThdCount = 32;
		if ((iEnd - iBegin) < iTotalThdCount)
			iTotalThdCount = iEnd - iBegin;

		uint64_t iCounter = (iEnd - iBegin) / iTotalThdCount;



		for (uint64_t i = 0; i != iTotalThdCount; ++i)
		{
			uint64_t iNext = iBegin + iCounter;
			thread thd{ &CPrimeNr::FindPrimeNumbers,this ,iBegin, iNext };
			vThd.push_back(move(thd));
			iBegin = iNext;
		}
		if (iBegin != iEnd)
		{
			thread thd{ &CPrimeNr::FindPrimeNumbers,this ,iBegin, iEnd };
			vThd.push_back(move(thd));
		}

		for (uint32_t i = 0; i != vThd.size(); ++i)
		{
			vThd[i].join();
		}

	}
};



class CPrimeNr_SimpleMethod : public CPrimeNr
{
public:
	bool IsPrime(_int_t iNr)
	{
		if (iNr <= 1)
			return false;

		if (iNr <= 3)
			return true;
		else if (!(iNr % 2) || !(iNr % 3))
			return false;

		

		_int_t iEnd = (_int_t)sqrt(iNr) + 1;

		for (_int_t i = 5;i <= iEnd;i += 6)
		{
			if (!(iNr%i) || !(iNr % (i + 2)))
				return false;
		}
		return true;
	}
};


class CPrimeNr_Fermat : public CPrimeNr
{
public:
	CPrimeNr_Fermat(void)
	{
		srand((unsigned int)time(0));
	}
	bool IsPrime(_int_t iNr)
	{
		return Fermat(iNr, 50);
	}

private:

	_int_t modulo(_int_t base, _int_t exponent, _int_t mod)
	{
		_int_t x = 1;
		_int_t y = base;
		while (exponent > 0)
		{
			if (exponent % 2 == 1)
				x = (x * y) % mod;
			y = (y * y) % mod;
			exponent = exponent / 2;
		}
		return x % mod;
	}

	bool Fermat(_int_t p, int iterations)
	{
		if (p == 1)
		{
			return false;
		}
		for (int i = 0; i < iterations; i++)
		{
			_int_t a = rand() % (p - 1) + 1;
			if (modulo(a, p - 1, p) != 1)
			{
				return false;
			}
		}
		return true;
	}
};

class CPrimeNr_Miller : public CPrimeNr
{
public:
	CPrimeNr_Miller(void)
	{
		srand((unsigned int)time(0));
	}
	
	_int_t modulo(_int_t base, _int_t exponent, _int_t mod)
	{
		_int_t x = 1;
		_int_t y = base;
		while (exponent > 0)
		{
			if (exponent % 2 == 1)
				x = (x * y) % mod;
			y = (y * y) % mod;
			exponent = exponent / 2;
		}
		return x % mod;
	}
	_int_t mulmod(_int_t a, _int_t b, _int_t mod)
	{
		_int_t x = 0, y = a % mod;
		while (b > 0)
		{
			if (b % 2 == 1)
			{
				x = (x + y) % mod;
			}
			y = (y * 2) % mod;
			b /= 2;
		}
		return x % mod;
	}

	bool IsPrime(_int_t n)
	{
		if (n <= 1)
			return false;

		if (n <= 3)
			return true;

		if (!(n % 2))
			return false;

		register _int_t d = n - 1;
		register _int_t r = 0;

		while (!(d&1))
		{
			d >>= 1;
			r++;
		}


		register int k = 50;

		for (register int i = 0; i != k; ++i)
		{
			_int_t a = rand() % (n - 3) + 2;
			//_int_t a = rand() % (n - 1) + 1;
			//_int_t x = Pow(a, d) % n;
			_int_t x = modulo(a, d, n);
/*
			if (x1 != x)
			{
				int dupa = 55;
				dupa++;
			}*/

			if ((x == 1) || (x == (n-1)))
				continue;

			if (InLoop(x, n,r))
				continue;

			return false;
		}
		return true;
	}

private:
	bool InLoop(_int_t &x, _int_t &n, register _int_t r)
	{
		while (--r)
		{
			_int_t x1 = mulmod(x, x, n);

			x = (x*x) % n;

			if (x1 != x)
			{
				int dupa = 55;
				dupa++;
			}


			if (x == 1)
				return false;

			if (x == n-1)
				return true;
		}
		return false;
	}
	_int_t Pow(_int_t x, _int_t n)
	{

		if (!n)
			return 1;

		if (n&1)
		{
			_int_t iTmp = Pow(x, ((n - 1)/2));
			return x*iTmp*iTmp;
		}
		else
		{
			_int_t iTmp = Pow(x, (n / 2));
			return iTmp*iTmp;
		}
	}


};

class CPrimeNr_Rabin : public CPrimeNr
{
public:
	CPrimeNr_Rabin(void)
	{
		srand((unsigned int)time(0));
	}

	using ll = _int_t;

	bool IsPrime(_int_t iNr)
	{
		return Miller(iNr, 50);
	}


private:
	ll mulmod(ll a, ll b, ll mod)
	{
		ll x = 0, y = a % mod;
		while (b > 0)
		{
			if (b % 2 == 1)
			{
				x = (x + y) % mod;
			}
			y = (y * 2) % mod;
			b /= 2;
		}
		return x % mod;
	}
	/*
	 * modular exponentiation
	  */
	ll modulo(ll base, ll exponent, ll mod)
	{
		ll x = 1;
		ll y = base;
		while (exponent > 0)
		{
			if (exponent % 2 == 1)
				x = (x * y) % mod;
			y = (y * y) % mod;
			exponent = exponent / 2;
		}
		return x % mod;
	}

	/*
	 * Miller-Rabin primality test, iteration signifies the accuracy
	  */
	bool Miller(ll p, int iteration)
	{
		if (p < 2)
		{
			return false;
		}
		if (p != 2 && p % 2 == 0)
		{
			return false;
		}
		ll s = p - 1;
		while (s % 2 == 0)
		{
			s /= 2;
		}
		for (int i = 0; i < iteration; i++)
		{
			ll a = rand() % (p - 1) + 1, temp = s;
			ll mod = modulo(a, temp, p);
			while (temp != p - 1 && mod != 1 && mod != p - 1)
			{
				mod = mulmod(mod, mod, p);
				temp *= 2;
			}
			if (mod != p - 1 && temp % 2 == 0)
			{
				return false;
			}
		}
		return true;
	}
};


bool IsPrimaryNumber(uint64_t iNr)
{
	if (iNr < 2)
		throw string("Wrong Nr");

	if (!(iNr % 2))
		return false;

	uint64_t iEnd = (uint64_t)sqrt((double)iNr)+1;
	for (uint64_t i = 3; i <= iEnd; i+=2)
	{
		if (!(iNr%i))
			return false;
	}
	return true;
}

void OnPrimaryNumberFound(uint64_t iNr)
{
	static mutex mMutex;
	static int iCounter = 0;

	mMutex.lock();
	cout << "Liczba Pierwsza[" << iCounter++ <<"]: "<< iNr << endl;
	mMutex.unlock();

}

void FindPrimaryNumbers(uint64_t iBegin, uint64_t iEnd)
{
	if (iBegin < 2)
		iBegin = 3;
	if (iBegin > iEnd)
		return;

	if (!(iBegin % 2))
		iBegin++;

	for (uint64_t i = iBegin;i <= iEnd; i+=2)
	{
		if (IsPrimaryNumber(i))
		{
			OnPrimaryNumberFound(i);
		}
	}
}


void MultiThread_FindPrimaryNumbers(uint64_t iBegin, uint64_t iEnd)
{
	vector<thread> vThd;
	uint64_t iTotalThdCount = 16;
	if ((iEnd - iBegin) < iTotalThdCount)
		iTotalThdCount = iEnd - iBegin;

	uint64_t iCounter = (iEnd - iBegin) / iTotalThdCount;
	for (uint64_t i = 0;i != iTotalThdCount; ++i)
	{
		uint64_t iNext = iBegin + iCounter;
		thread thd{ FindPrimaryNumbers,iBegin, iNext };
		vThd.push_back(move(thd));
		iBegin = iNext;
	}
	if (iBegin != iEnd)
	{
		thread thd{ FindPrimaryNumbers,iBegin, iEnd };
		vThd.push_back(move(thd));
	}

	for (uint32_t i = 0; i != vThd.size(); ++i)
	{
		vThd[i].join();
	}


}


class CPrime_Strong
{
public:
	bool IsPrime(_int_t iNr)
	{
		CPrimeNr_SimpleMethod mMth;

		if (mMth.IsPrime(iNr))
		{
			iNr = (iNr-1)/2;
			
			return mMth.IsPrime(iNr);

		}
		return false;
	}
};

class CFactoring_Simple
{
public:

	vector<_int_t> Factor(_int_t iNr)
	{
		vector<_int_t> viTmp;
		_int_t iResp = FactorBase(iNr);
		viTmp.push_back(iResp);
		if (iResp == iNr)
		{
			return viTmp;
		}
		else
		{
			auto viResp = Factor(iNr / iResp);

			for (int i = 0; i != viResp.size(); ++i)
			{
				viTmp.push_back(viResp[i]);
			}
		}

		return viTmp;
	}




private:
	_int_t FactorBase(_int_t iNr)
	{
		_int_t iEnd = (_int_t)sqrt(iNr) + 1;
		while (iEnd*iEnd < iNr)
			iEnd++;

		for (_int_t i = 5; i <= iEnd; i += 6)
		{
			if (!(iNr%i))
				return i;
				
			if (!(iNr % (i + 2)))
			{
				return i + 2;
				break;
			}
		}

		return iNr;
	}
};



void Factoring(_int_t iVal)
{
	

	CFactoring_Simple mFactor;

	auto t1 = steady_clock::now();
	auto viFactor = mFactor.Factor(iVal);



	auto d = steady_clock::now() - t1;
	auto ms = duration_cast<milliseconds>(d).count();

	cout << "Czas: " << (double)ms / 1000 << endl;


	for (uint32_t i = 0; i != viFactor.size(); ++i)
	{
		cout << "LF[" << i << "]: " << viFactor[i] << endl;
	}
}


void Factoring(void)
{
	CPrime_Strong mSimple;
	vector<_int_t> m_viPrime;
	for (_int_t i = 100000000; ; ++i)
	{
		if (mSimple.IsPrime(i))
		{
			m_viPrime.push_back(i);
			if (m_viPrime.size() == 2)
				break;
		}
	}

	cout << "Wygenerowane liczby" << endl;
	for (uint32_t i = 0; i != m_viPrime.size(); ++i)
	{
		cout << "L[" << i << "]: " << m_viPrime[i] << endl;
	}
	_int_t iMulNr = m_viPrime[0] * m_viPrime[1];
	//iMulNr += 2;
	cout << "Iloczyn liczb: " << iMulNr << endl;


	Factoring(iMulNr);
}


void FindPrime(void)
{
	cout << "Wybierz algorytm wyszukiwania" << endl;


	string mStr;
	bool bFirst = true;

	do
	{

		if (!bFirst)
		{
			system("cls");
			cout << "Wprowadzono niewlasciwa wartosc, powtorz" << endl;
		}


		bFirst = false;

		cout << "1: Metoda Naiwna" << endl;
		cout << "2: Algorytm Millera" << endl;
		cout << "3: Algorytm Fermata" << endl;
		cout << "E - konczy prace programu" << endl;

		char pTmp[1024];
		memset(pTmp, 0, sizeof(pTmp));

		cin >> pTmp;

		mStr = string(pTmp);




	} while (mStr.length() != 1 || !(mStr[0] == '1' || mStr[0] == '2' || mStr[0] == '3' || mStr[0] == 'E' || mStr[0] == 'e'));

	CPrimeNr *pPrime = nullptr;
	switch (mStr[0])
	{
	case '1':
		pPrime = new CPrimeNr_SimpleMethod();
		break;
	case '2':
		pPrime = new CPrimeNr_Miller();
		break;
	case '3':
		pPrime = new CPrimeNr_Fermat();
		break;
	case 'E':
	case 'e':
		exit(0);
		break;
	}


	

	_int_t iStart = 2305843009213693951;
	_int_t iEnd = iStart + 5000;

/*
	cout << "Wprowadz wartosc poczatkowa wyszukiwania" << endl;
	cin >> iStart;

	cout << "Wprowadz wartosc koncowa wyszukiwania" << endl;
	cin >> iEnd;
	*/


	

	pPrime->Multithread_FindPrimaryNumbers(iStart, iEnd);

	

	delete pPrime;
}

//labelFor:
typedef int64_t _uint_t;
vector<_uint_t> ResztaKwadratowa(_uint_t iP)
{
	vector<_uint_t> mTmp;
	for (_uint_t i = 1; i != iP; ++i)
	{
		
		_uint_t iNr = i*i;
		iNr %= iP;
		cout << "Reszta kwadratowa [" << i << "] = " << iNr << endl;
		for (_uint_t iCnt = 0; iCnt != mTmp.size(); iCnt++)
		{

			if (iNr == mTmp[iCnt])
			{
				goto labelFor;
			}
		}
		
		mTmp.push_back(iNr);


	labelFor:;

	}

	return mTmp;
}

uint32_t FastPow(_uint_t iNr, _uint_t iPow)
{
	_uint_t iResp = 1;
	for (_uint_t i = 0; i < iPow; ++i)
	{
		iResp *= iNr;
	}

	return iResp;
}



_uint_t PowMod(_uint_t base, _uint_t exponent, _uint_t mod)
{
	_uint_t x = 1;
	_uint_t y = base;
	while (exponent > 0)
	{
		if (exponent % 2 == 1)
			x = (x * y) % mod;
		y = (y * y) % mod;
		exponent = exponent / 2;
	}
	return x % mod;
}

_uint_t Legendre(_uint_t a, _uint_t n)
{
	_uint_t iRes = PowMod(a, (n - 1) / 2, n);


	if (iRes == n - 1)
	{
		return -1;
	}

	return iRes;
}


void Check_Legendre(void)
{
	//_uint_t Legendre(_uint_t a, _uint_t n)
	_uint_t n = 519;

/*
	for (_uint_t i = 510; i < 1000; ++i)
	{
		if (IsPrimaryNumber(i))
		{
			n = i;
			break;
		}
	}
*/
	bool bIsPrime = IsPrimaryNumber(n);

	if (bIsPrime)
	{
		cout << "Liczba jest pierwsza " << n << endl;
	}
	else
	{
		cout << "Liczba nie jest pierwsza " << n << endl;
	}

	for (_uint_t i = 3; i < n; i += 1)
	{
		_uint_t iLegendre = Legendre(i, n);


		cout << "Legendre [" << i << "]: " << iLegendre << endl;

	}


}


void ResztyKwadratowe(void)
{
	Check_Legendre();
	return;
	_uint_t iNr = 17;
	vector<_uint_t> mData = ResztaKwadratowa(iNr);
	
	cout << "Reszty kwadratowe: dla: p="<< iNr << endl;

	for (_uint_t i = 0; i != mData.size(); ++i)
	{
		cout << "Reszty kwadratowe: "<< mData[i] << endl;
	}

	for (_uint_t i = 1; i != iNr; ++i)
	{
		_uint_t iRes = FastPow(i, (iNr - 1) / 2) % iNr;

		iRes = PowMod(i, (iNr - 1) / 2, iNr);
		cout << "Wynik Operacji dla :" << i << ": " << iRes << endl;
	}

}




typedef struct
{
	uint64_t m_iP;
	uint64_t m_iQ;
	uint64_t m_iE;
	uint64_t m_iD;
	uint64_t m_iN;
}rsa_t;

uint64_t Nwd(uint64_t a, uint64_t b)
{
	uint64_t c;                    // Tworzymy zmienną c o typie int
	while (b != 0)             // Tworzymy pętlę działającą gdy b ≠ 0
	{
		c = a % b;              // Zmienna c przyjmuje wartość a modulo b
		a = b;                // Przypisujemy b do a
		b = c;                // Przypisujemy c do b
	}
	return a;                 // Zwracamy a, czyli żądane NWD(a,b)

}

uint64_t Extended_Euklides(uint64_t nwd_b, uint64_t nwd_a)
{
	uint64_t r, nwd, a, q, b;
	uint64_t x, x1, x2;
	uint64_t y, y1, y2;




	// a must be greater than b
	if (nwd_b > nwd_a)
	{
		nwd = nwd_b;
		nwd_b = nwd_a;
		nwd_a = nwd;
	}


	//initialize a and b
	a = nwd_a;
	b = nwd_b;

	//initialize r and nwd
	q = a / b;
	r = a - q*b;
	nwd = b;

	//initialize x and y
	x2 = 1;
	x1 = 0;
	y2 = 0;
	y1 = 1;
	x = 1;
	y = y2 - (q - 1)*y1;

	while (r != 0)
	{
		a = b;
		b = r;

		x = x2 - q*x1;
		x2 = x1;
		x1 = x;

		y = y2 - q*y1;
		y2 = y1;
		y1 = y;

		nwd = r;
		q = a / b;
		r = a - q*b;
	}
	return y;
}



void Rsa_Init(rsa_t &rRsa)
{
	srand(time(0));
	uint32_t iNr;
	uint64_t iFi;
	do
	{
		memset(&rRsa, 0, sizeof(rsa_t));

		while (1)
		{
			iNr = (uint16_t)rand();
			//iNr |= rand()<<16;
			if (IsPrimaryNumber(iNr))
			{
				rRsa.m_iP = iNr;
				break;
			}
		}
		while (1)
		{
			iNr = (uint16_t)rand();
			//iNr |= rand() << 16;
			if (IsPrimaryNumber(iNr))
			{
				rRsa.m_iQ = iNr;
				break;
			}
		}

		rRsa.m_iN = rRsa.m_iQ * rRsa.m_iP;


		iFi = (rRsa.m_iQ - 1)*(rRsa.m_iP - 1);
		for (uint64_t i = rRsa.m_iQ + 1; i < iFi; ++i)
		{
			uint64_t iTmp = Nwd(i, iFi);

			if (iTmp == 1)
			{
				rRsa.m_iE = i;
				break;
			}

		}

		uint64_t iOdw = Extended_Euklides(rRsa.m_iE, iFi);

		rRsa.m_iD = iOdw;

		if ((rRsa.m_iD*rRsa.m_iE) % iFi == 1)
		{
			break;
		}


	} while (1);
	//cout << "Rsa Fi " << iFi << endl;
	cout << "Rsa E*DModFi " << (rRsa.m_iD*rRsa.m_iE)% iFi << endl;
}




/*
void Faktoryzacja(void)
{
	uint64_t iNr = 0;
	cout << "Podaj liczbę do faktoryzacji" << endl;
	cin >> iNr;

	CFactoring_Simple mFactor;

	auto t1 = steady_clock::now();
	auto viFactor = mFactor.Factor(iVal);



	auto d = steady_clock::now() - t1;
	auto ms = duration_cast<milliseconds>(d).count();

	cout << "Czas: " << (double)ms / 1000 << endl;


	for (uint32_t i = 0; i != viFactor.size(); ++i)
	{
		cout << "LF[" << i << "]: " << viFactor[i] << endl;
	}



}
*/


void Rsa_EnCrypt(vector<uint64_t> &rData, rsa_t &rRsa)
{
	for (auto &p : rData)
	{
		p = PowMod(p, rRsa.m_iE, rRsa.m_iN);
	}
}
void Rsa_DeCrypt(vector<uint64_t> &rData, rsa_t &rRsa)
{
	for (auto &p : rData)
	{
		p = PowMod(p, rRsa.m_iD, rRsa.m_iN);
	}
}

void RSA_Test(void)
{
	rsa_t mRsa;
	Rsa_Init(mRsa);

	cout << "Rsa P " << mRsa.m_iP << endl;
	cout << "Rsa Q " << mRsa.m_iQ << endl;
	cout << "Rsa N " << mRsa.m_iN << endl;
	cout << "Rsa E " << mRsa.m_iE << endl;
	cout << "Rsa D " << mRsa.m_iD << endl;


	string mStr("Ala ma kota i go bardzo lubi");

	vector<uint64_t> rsaMsg;

	for (uint32_t i = 0; i != mStr.length(); ++i)
	{
		rsaMsg.push_back(mStr[i]);
		cout << "Msg D[" << i << "]: " << rsaMsg[i] << endl;
	}

	Rsa_EnCrypt(rsaMsg, mRsa);

	cout << endl << endl;
	for (uint32_t i = 0; i != rsaMsg.size(); ++i)
	{
		cout << "Msg C[" << i << "]: " << rsaMsg[i] << endl;
	}


	cout << endl << endl;
	Rsa_DeCrypt(rsaMsg, mRsa);
	for (uint32_t i = 0; i != rsaMsg.size(); ++i)
	{
		cout << "Msg DD[" << i << "]: " << rsaMsg[i] << endl;
	}

}

int main1()
{
	int r, nwd, a, q, b;
	int x, x1, x2;
	int y, y1, y2;
	int nwd_a, nwd_b;

	//get all data
	printf("Podaj pierwsza liczbe\n");
	scanf("%d", &nwd_a);

	printf("Podaj druga liczbe\n");
	scanf("%d", &nwd_b);

	// a must be greater than b
	if (nwd_b > nwd_a)
	{
		nwd = nwd_b;
		nwd_b = nwd_a;
		nwd_a = nwd;
	}


	//initialize a and b
	a = nwd_a;
	b = nwd_b;

	//initialize r and nwd
	q = a / b;
	r = a - q*b;
	nwd = b;

	//initialize x and y
	x2 = 1;
	x1 = 0;
	y2 = 0;
	y1 = 1;
	x = 1;
	y = y2 - (q - 1)*y1;

	while (r != 0)
	{
		a = b;
		b = r;

		x = x2 - q*x1;
		x2 = x1;
		x1 = x;

		y = y2 - q*y1;
		y2 = y1;
		y1 = y;

		nwd = r;
		q = a / b;
		r = a - q*b;
	}

	//present results
	printf("NWD(%d, %d) = %d = %d * %d + %d * %d\n", nwd_a, nwd_b, nwd, x, nwd_a, y, nwd_b);

	if (nwd == 1)
		printf("%d * %d mod %d = 1", nwd_b, y, nwd_a);
	system("pause");
	return 0;
}
int main(void)
{
	RSA_Test();
	//ResztyKwadratowe();

	system("pause");
	return 0;

	Factoring();
	system("pause");
	return 0;

	FindPrime();
	system("pause");
	return 0;
	cout << "Czesc" << endl;
	uint64_t iBegin = 2305843009213693951;
	uint64_t iEnd = iBegin + 160;

	CPrimeNr_SimpleMethod mSimple;
	CPrimeNr_Miller mMiller;
	CPrimeNr_Rabin mRabin;

	


	for (_int_t i = 2; i != 10000; ++i)
	{
		bool b1 = mSimple.IsPrime(i);
		bool b2 = mMiller.IsPrime(i);
		bool b3 = mRabin.IsPrime(i);
		if (b1 != b2 || b3 != b2)
		{
			cout << "Rozbierzność: " << i << "\t\t B1:" << b1 << " \tB2: " << b2 << " \tB3: " << b3 << endl;
		}


	}
	system("pause");
	return 0;

	/*
	for (int x = 1; x != 20; ++x)
	{
		for (int y = 0; y != 10; ++y)
		{
			_int_t iRes1 = Pow(x,y);
			_int_t iRes2 = pow(x, y);

			cout <<"x:"<<x<<"y: "<<y <<"\t\t" << "Res1: " << iRes1 << "\t\t" << "Res2: " << iRes2 << endl;

			if (iRes1 != iRes2)
			{
				system("pause");
				return 0;
			}

		}
	}
	system("pause");
	return 0;
*/
	cout << "Fast search" << endl;
	
	steady_clock::time_point t1;
	int64_t ms;
	steady_clock::duration d;
	t1 = steady_clock::now();
	MultiThread_FindPrimaryNumbers(iBegin, iEnd);
	d = steady_clock::now() - t1; 	ms = duration_cast<milliseconds>(d).count();
	cout << "TimeSpend[s]: " << (double)ms / 1000.0f << endl;
/*
	cout << "Slow search" << endl;
	t1 = steady_clock::now();
	FindPrimaryNumbers(iBegin, iEnd);
	d = steady_clock::now() - t1; 	ms = duration_cast<milliseconds>(d).count();
	cout << "TimeSpend[s]: " << (double)ms / 1000.0f << endl;
	*/

/*
	cout << "Class search" << endl;
	t1 = steady_clock::now();
	CPrimeNr mPrimeNr;
	mPrimeNr.Multithread_FindPrimaryNumbers(iBegin, iEnd);
	d = steady_clock::now() - t1; 	ms = duration_cast<milliseconds>(d).count();
	cout << "TimeSpend[s]: " << (double)ms / 1000.0f << endl;

	*/
	system("pause");
	return 0;
}











