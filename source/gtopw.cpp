#include "../headers/gtopw.h"
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include "../headers/auxfun1.h"
#include "../headers/auxfun2.h"
#include "../headers/caps.h"
#include "omp.h"

using namespace std;

int belt[500];

double fact[161];
double dfact[201];
double binom[101][101];
double omega[33][33];

double clmr[bas_lmax + 1][2 * bas_lmax + 1][crt_siz[bas_lmax]] = {{{0.0}}};
double xyz_norm[bas_lmax + 1][crt_siz[bas_lmax]] = {{0.0}};

<<<<<<< HEAD
int main(int argc, char *argv[]) {
	/* put a battery in the clock */
	double t_start, t_end;
	t_start = omp_get_wtime();

	/* code initialisation */
	Keywords keys;
	Init();
	string inpname;
	if (argc > 1) {
		inpname = argv[1];
	} else {
		cout << " Name of the input file not specified! Emergency halt." << endl;
		exit(EXIT_FAILURE);
	};

	/* read the $BASIS keyword */
	vector<GaussC> basis;
	vector<Nuclei> nucli;

	GaussC gau_tmp;
	Nuclei nuc_tmp;

	int num_crt, num_sph;
	num_crt = 0;
	num_sph = 0;

	ifstream ifile;
	string line, name;
	stringstream ss;

	/* read destination path (MS)*/
	ifile.open(inpname);
	while (getline(ifile, line))
		if (line == "$PATH") {
			getline(ifile, line);
			keys.path = line;
			break;
		};
	ifile.close();
	line.clear();

	keys.name_me(inpname);

	ifile.open(inpname);

	do {
		getline(ifile, line);
	} while (line != "$BASIS");
	getline(ifile, line);
	do {
		ss.clear();
		ss << line;
		line.clear();
		ss >> nuc_tmp.name >> nuc_tmp.chrg >> nuc_tmp.Cx >> nuc_tmp.Cy >> nuc_tmp.Cz;

		nucli.push_back(std::move(nuc_tmp));
		getline(ifile, line);
		string dummy, moment;
		int end;

		double ex, d_re, d_im;
		do {
			ss.clear();
			ss << line;
			line.clear();
			ss >> moment >> end;

			if (moment == "S")
				gau_tmp.lA = 0;
			else if (moment == "P")
				gau_tmp.lA = 1;
			else if (moment == "D")
				gau_tmp.lA = 2;
			else if (moment == "F")
				gau_tmp.lA = 3;
			else if (moment == "G")
				gau_tmp.lA = 4;
			else if (moment == "H")
				gau_tmp.lA = 5;
			else if (moment == "I")
				gau_tmp.lA = 6;
			else if (moment == "K")
				gau_tmp.lA = 7;
			else if (moment == "L")
				gau_tmp.lA = 8;
			else {
				cout << " Unrecognised angular momentum of a function!" << endl;
				cout << " Emergency halt." << endl;
				exit(EXIT_FAILURE);
			};

			num_crt += crt_siz[gau_tmp.lA];
			num_sph += sph_siz[gau_tmp.lA];

			gau_tmp.clen = end;
			for (int i = 0; i < end; i++) {
				getline(ifile, line);
				ss.clear();
				ss << line;
				line.clear();
				ss >> dummy >> ex >> d_re >> d_im >> gau_tmp.kx >> gau_tmp.ky >> gau_tmp.kz;
				gau_tmp.alphaA.push_back(ex);
				gau_tmp.dA_re.push_back(d_re);
				gau_tmp.dA_im.push_back(d_im);
			};

			RenormContr(gau_tmp, keys.norm1E);
			gau_tmp.Ax = nuc_tmp.Cx;
			gau_tmp.Ay = nuc_tmp.Cy;
			gau_tmp.Az = nuc_tmp.Cz;

			basis.push_back(std::move(gau_tmp));
			getline(ifile, line);
		} while (!line.empty());
		getline(ifile, line);
	} while (line != "$END");
	ifile.close();

	/* basis set data - mostly optional printing */
	cout << endl;
	cout << " Number of basis set shells      = " << basis.size() << endl;
	cout << " Number of functions (spherical) = " << num_sph << endl;
	cout << " Number of functions (cartesian) = " << num_crt << endl;
	cout << " Contracted functions are now normalised to the unity." << endl;
	cout << endl;

	/* read the remaining keywords */
	ifile.open(inpname);
	while (getline(ifile, line))
		if (line == "$INTS")
			break;
	getline(ifile, line);
	int num_ints = atoi(line.c_str());

	for (int i = 0; i < num_ints; i++) {
		getline(ifile, line);
		if (line == "STVH")
			keys.stvh = true;
		else if (line == "DIPOLE")
			keys.dip = true;
		else if (line == "QUADRUP")
			keys.quad = true;
		else if (line == "VELOCITY")
			keys.velo = true;
		else if (line == "ERI")
			keys.eri = true;
		else if (line == "PROJ")
			keys.proj = true;
		else if (line == "KINXYZ")
			keys.kinc = true;
		else if (line == "CAPINT")
			keys.capi = true;
		else
			cout << " Unrecognised integral type (" << i + 1 << "). Omitted." << endl;
	};

	cout << " Requested one-electron integrals:" << endl;
	if (keys.stvh)
		cout << " - ordinary one-electron STVH integrals" << endl;
	if (keys.dip)
		cout << " - dipole moment integrals" << endl;
	if (keys.quad)
		cout << " - quadruple moment integrals" << endl;
	if (keys.velo)
		cout << " - dipole velocity integrals" << endl;
	if (keys.proj)
		cout << " - projection integrals" << endl;
	if (keys.kinc)
		cout << " - kinetic energy components int." << endl;
	if (keys.capi)
		cout << " - complex absorbing potentials int." << endl;
	cout << endl;
	ifile.close();

	int num_points(0);
	line.clear();
	ifile.open(inpname);
	while (getline(ifile, line))
		if (line == "$POINTS") {
			getline(ifile, line);
			num_points = atoi(line.c_str());
			break;
		};

	cout << " Points for properties calculations (" << num_points << "): " << endl;
	for (int i = 0; i < num_points; i++) {
		getline(ifile, line);
		ss.clear();
		ss << line;
		double xp, yp, zp;
		ss >> xp >> yp >> zp;
		keys.points.push_back(std::make_tuple(xp, yp, zp));
		cout << " " << std::get<0>(keys.points[i]) << " " << std::get<1>(keys.points[i]) << " " << std::get<2>(keys.points[i]) << endl;
	};
	ifile.close();

	line.clear();
	ifile.open(inpname);
	while (getline(ifile, line))
		if (line == "$RGRID") {
			getline(ifile, line);
			keys.n_rad = atoi(line.c_str());
			getline(ifile, line);
			keys.rad_str = atof(line.c_str());
			getline(ifile, line);
			keys.rad_stp = atof(line.c_str());
			getline(ifile, line);
			keys.lproj_max = atoi(line.c_str());
			break;
		};

	cout << endl;
	if (keys.n_rad <= 0 ||
	    keys.rad_stp <= 0.0 ||
	    keys.rad_str < 0.0)
		keys.proj = false;  /// fail-safe

	if (keys.proj) {
		cout << " Radial grid parameters:" << endl;
		cout << "  Number of radial points     = " << keys.n_rad << endl;
		cout << "  Maximal projection momentum = " << keys.lproj_max << endl;
		cout << "  Starting point of the grid  = " << keys.rad_str << endl;
		cout << "  Grid step parameters        = " << keys.rad_stp << endl;
	};
	cout << endl;
	ifile.close();

	line.clear();
	ifile.open(inpname);
	while (getline(ifile, line))
		if (line == "$CAPAR") {
			getline(ifile, line);
			keys.cap_on = atof(line.c_str());
			getline(ifile, line);
			keys.cap_mul = atof(line.c_str());
			break;
		};

	if (keys.cap_on < 0.0 ||
	    keys.cap_mul < 0.0)
		keys.capi = false;  /// fail-safe

	if (keys.capi) {
		cout << " Complex absorbing potentials parameters:" << endl;
		cout << "  Onset of the potential   = " << keys.cap_on << endl;
		cout << "  Multiplicative parameter = " << keys.cap_mul << endl;
	};
	cout << endl;
	ifile.close();

	/* ---------------------------------- */
	/* ----- ONE-ELECTRON INTEGRALS ----- */
	/* ---------------------------------- */

	int num_crt2, num_sph2;
	int shl_siz_crt, shl_siz_crt2, shl_siz_sph, shl_siz_sph2;
	num_crt2 = num_crt * num_crt;
	num_sph2 = num_sph * num_sph;

	shl_siz_crt = crt_siz[bas_lmax];
	shl_siz_sph = sph_siz[bas_lmax];
	shl_siz_crt2 = shl_siz_crt * shl_siz_crt;
	shl_siz_sph2 = shl_siz_sph * shl_siz_sph;

	/* memory for complete matrices */
	arr1d<cdouble> ovrl_crt(num_crt2);
	arr1d<cdouble> ovrl_sph(num_crt2);

	arr1d<cdouble> kin_crt(num_crt2);
	arr1d<cdouble> kin_sph(num_crt2);

	arr1d<cdouble> nuc_crt(num_crt2);
	arr1d<cdouble> nuc_sph(num_crt2);

	arr1d<cdouble> bare_crt(num_crt2);
	arr1d<cdouble> bare_sph(num_crt2);

	arr1d<cdouble> dipx_crt(num_crt2);
	arr1d<cdouble> dipx_sph(num_crt2);

	arr1d<cdouble> dipy_crt(num_crt2);
	arr1d<cdouble> dipy_sph(num_crt2);

	arr1d<cdouble> dipz_crt(num_crt2);
	arr1d<cdouble> dipz_sph(num_crt2);

	arr1d<cdouble> qdxx_crt(num_crt2);
	arr1d<cdouble> qdxx_sph(num_crt2);

	arr1d<cdouble> qdyy_crt(num_crt2);
	arr1d<cdouble> qdyy_sph(num_crt2);

	arr1d<cdouble> qdzz_crt(num_crt2);
	arr1d<cdouble> qdzz_sph(num_crt2);

	arr1d<cdouble> qdxy_crt(num_crt2);
	arr1d<cdouble> qdxy_sph(num_crt2);

	arr1d<cdouble> qdxz_crt(num_crt2);
	arr1d<cdouble> qdxz_sph(num_crt2);

	arr1d<cdouble> qdyz_crt(num_crt2);
	arr1d<cdouble> qdyz_sph(num_crt2);

	arr1d<cdouble> grdx_crt(num_crt2);
	arr1d<cdouble> grdx_sph(num_crt2);

	arr1d<cdouble> grdy_crt(num_crt2);
	arr1d<cdouble> grdy_sph(num_crt2);

	arr1d<cdouble> grdz_crt(num_crt2);
	arr1d<cdouble> grdz_sph(num_crt2);

	arr1d<cdouble> tx_crt(num_crt2);
	arr1d<cdouble> tx_sph(num_crt2);

	arr1d<cdouble> ty_crt(num_crt2);
	arr1d<cdouble> ty_sph(num_crt2);

	arr1d<cdouble> tz_crt(num_crt2);
	arr1d<cdouble> tz_sph(num_crt2);

	arr1d<cdouble> angz_crt(num_crt2);
	arr1d<cdouble> angz_sph(num_crt2);

	arr1d<cdouble> capi_crt(num_crt2);
	arr1d<cdouble> capi_sph(num_crt2);

	/* memory for shell data */
	arr1d<cdouble> shl_crt_s(shl_siz_crt2);
	arr1d<cdouble> shl_sph_s(shl_siz_sph2);

	arr1d<cdouble> shl_crt_t(shl_siz_crt2);
	arr1d<cdouble> shl_sph_t(shl_siz_sph2);

	arr1d<cdouble> shl_crt_v(shl_siz_crt2);
	arr1d<cdouble> shl_sph_v(shl_siz_sph2);

	arr1d<cdouble> shl_crt_dx(shl_siz_crt2);
	arr1d<cdouble> shl_sph_dx(shl_siz_sph2);

	arr1d<cdouble> shl_crt_dy(shl_siz_crt2);
	arr1d<cdouble> shl_sph_dy(shl_siz_sph2);

	arr1d<cdouble> shl_crt_dz(shl_siz_crt2);
	arr1d<cdouble> shl_sph_dz(shl_siz_sph2);

	arr1d<cdouble> shl_crt_xx(shl_siz_crt2);
	arr1d<cdouble> shl_sph_xx(shl_siz_sph2);

	arr1d<cdouble> shl_crt_yy(shl_siz_crt2);
	arr1d<cdouble> shl_sph_yy(shl_siz_sph2);

	arr1d<cdouble> shl_crt_zz(shl_siz_crt2);
	arr1d<cdouble> shl_sph_zz(shl_siz_sph2);

	arr1d<cdouble> shl_crt_xy(shl_siz_crt2);
	arr1d<cdouble> shl_sph_xy(shl_siz_sph2);

	arr1d<cdouble> shl_crt_xz(shl_siz_crt2);
	arr1d<cdouble> shl_sph_xz(shl_siz_sph2);

	arr1d<cdouble> shl_crt_yz(shl_siz_crt2);
	arr1d<cdouble> shl_sph_yz(shl_siz_sph2);

	arr1d<cdouble> shl_crt_gx(shl_siz_crt2);
	arr1d<cdouble> shl_sph_gx(shl_siz_sph2);

	arr1d<cdouble> shl_crt_gy(shl_siz_crt2);
	arr1d<cdouble> shl_sph_gy(shl_siz_sph2);

	arr1d<cdouble> shl_crt_gz(shl_siz_crt2);
	arr1d<cdouble> shl_sph_gz(shl_siz_sph2);

	arr1d<cdouble> shl_crt_tx(shl_siz_crt2);
	arr1d<cdouble> shl_sph_tx(shl_siz_sph2);

	arr1d<cdouble> shl_crt_ty(shl_siz_crt2);
	arr1d<cdouble> shl_sph_ty(shl_siz_sph2);

	arr1d<cdouble> shl_crt_tz(shl_siz_crt2);
	arr1d<cdouble> shl_sph_tz(shl_siz_sph2);

	arr1d<cdouble> shl_crt_lz(shl_siz_crt2);
	arr1d<cdouble> shl_sph_lz(shl_siz_sph2);

	arr1d<cdouble> shl_crt_ca(shl_siz_crt2);
	arr1d<cdouble> shl_sph_ca(shl_siz_sph2);

	arr1d<cdouble> shl_crt_dum(shl_siz_crt2);
	arr1d<cdouble> shl_sph_dum(shl_siz_sph2);

	/* basic stuff for projections, CAPs, etc. */
	int ltot_max = 30;
	arr2d<double> alpha_ij(ltot_max + 1, ltot_max + 1);
	BuildAlphaIJ(ltot_max, alpha_ij.v);

	/* stuff for various property integrals */
	constexpr int dim_hash = 7;
	cdouble hash_fill[dim_hash][dim_hash] = {0.0};
	int hash[dim_hash][3] = {0};

	hash[0][0] = -1;
	hash[0][1] = +0;
	hash[0][2] = +0;

	hash[1][0] = +0;
	hash[1][1] = -1;
	hash[1][2] = +0;

	hash[2][0] = +0;
	hash[2][1] = +0;
	hash[2][2] = -1;

	hash[3][0] = +1;
	hash[3][1] = +0;
	hash[3][2] = +0;

	hash[4][0] = +0;
	hash[4][1] = +1;
	hash[4][2] = +0;

	hash[5][0] = +0;
	hash[5][1] = +0;
	hash[5][2] = +1;

	hash[6][0] = +0;
	hash[6][1] = +0;
	hash[6][2] = +0;

	int li, lj, lij, mi, mj, mjmx;
	int indi, indj;

	int posimx_crt, posjmx_crt;
	int posimx_sph, posjmx_sph;
	int shgi, shgj, shi, shj;
	int ia, ja, ka, ib, jb, kb;
	int ia1, ja1, ka1, ib1, jb1, kb1;

	double dum, ai, aj, aij;
	double Ax, Ay, Az, Bx, By, Bz;
	double ABx, ABy, ABz;
	double kPx, kPy, kPz, kP2;
	double Px, Py, Pz, kPP;
	double APx, APy, APz;
	double BPx, BPy, BPz;

	cdouble di, dj, dij, II;
	cdouble fac, Ex0, Ey0, Ez0;
	cdouble Ex0_1, Ey0_1, Ez0_1, Exyz;
	cdouble tmp1s, tmp2s, tmp3s, tmp4s;

	dum = 0.0;

	posimx_crt = 0;
	posimx_sph = 0;
	/* first loop over the shells */
	for (usint i = 0; i < basis.size(); i++) {
		li = basis[i].lA;
		Ax = basis[i].Ax;
		Ay = basis[i].Ay;
		Az = basis[i].Az;

		shi = sph_siz[li];
		shgi = crt_siz[li];

		posjmx_crt = 0;
		posjmx_sph = 0;
		/* second loop over the shells */
		for (usint j = 0; j <= i; j++) {
			lj = basis[j].lA;
			Bx = basis[j].Ax;
			By = basis[j].Ay;
			Bz = basis[j].Az;

			kPx = basis[j].kx - basis[i].kx;
			kPy = basis[j].ky - basis[i].ky;
			kPz = basis[j].kz - basis[i].kz;
			kP2 = kPx * kPx + kPy * kPy + kPz * kPz;

			ABx = Ax - Bx;
			ABy = Ay - By;
			ABz = Az - Bz;

			/* momentum transfer coefficients */
			Ex0 = exp(I * (basis[i].kx * Ax - basis[j].kx * Bx));
			Ey0 = exp(I * (basis[i].ky * Ay - basis[j].ky * By));
			Ez0 = exp(I * (basis[i].kz * Az - basis[j].kz * Bz));

			lij = li + lj;
			shj = sph_siz[lj];
			shgj = crt_siz[lj];

			ECoefs<double> Eijx(li + 1, lj + 1, dum, dum, ABx);
			ECoefs<double> Eijy(li + 1, lj + 1, dum, dum, ABy);
			ECoefs<double> Eijz(li + 1, lj + 1, dum, dum, ABz);

			RInts1E<cdouble> ROvrl(lij + 2);
			RInts1E<cdouble> RNucA(lij);

			/* CAP factors */
			arr2d<double> cap_ang(lij + 1, (lij + 1) * (lij + 2) / 2);
			arr2d<double> cap_rad(lij + 1, lij + 5);
			arr1d<cdouble> cap_tmp((lij + 1) * (lij + 2) / 2);

			if (keys.capi)
				CAPAngFac(kPx, kPy, kPz, lij, lij, alpha_ij.v, cap_ang.v);

			/* zero out the shell data */
			shl_crt_s.zero();
			shl_sph_s.zero();
			shl_crt_t.zero();
			shl_sph_t.zero();
			shl_crt_v.zero();
			shl_sph_v.zero();

			shl_crt_dx.zero();
			shl_sph_dx.zero();
			shl_crt_dy.zero();
			shl_sph_dy.zero();
			shl_crt_dz.zero();
			shl_sph_dz.zero();

			shl_crt_xx.zero();
			shl_sph_xx.zero();
			shl_crt_yy.zero();
			shl_sph_yy.zero();
			shl_crt_zz.zero();
			shl_sph_zz.zero();
			shl_crt_xy.zero();
			shl_sph_xy.zero();
			shl_crt_xz.zero();
			shl_sph_xz.zero();
			shl_crt_yz.zero();
			shl_sph_yz.zero();

			shl_crt_gx.zero();
			shl_sph_gx.zero();
			shl_crt_gy.zero();
			shl_sph_gy.zero();
			shl_crt_gz.zero();
			shl_sph_gz.zero();

			shl_crt_tx.zero();
			shl_sph_tx.zero();
			shl_crt_ty.zero();
			shl_sph_ty.zero();
			shl_crt_tz.zero();
			shl_sph_tz.zero();

			shl_crt_lz.zero();
			shl_sph_lz.zero();
			shl_crt_ca.zero();
			shl_sph_ca.zero();

			/* loop over contractions - shell I */
			for (usint k1 = 0; k1 < basis[i].alphaA.size(); k1++) {
				ai = basis[i].alphaA[k1];
				di = basis[i].dA_re[k1] - I * basis[i].dA_im[k1];

				/* loop over contractions - shell J */
				for (usint k2 = 0; k2 < basis[j].alphaA.size(); k2++) {
					aj = basis[j].alphaA[k2];
					dj = basis[j].dA_re[k2] + I * basis[j].dA_im[k2];

					aij = ai + aj;
					dij = di * dj;

					Px = ai * Ax + aj * Bx;
					Px = Px / aij;
					Py = ai * Ay + aj * By;
					Py = Py / aij;
					Pz = ai * Az + aj * Bz;
					Pz = Pz / aij;

					Eijx.zero();
					Eijy.zero();
					Eijz.zero();

					Eijx.load(ai, aj);
					Eijy.load(ai, aj);
					Eijz.load(ai, aj);

					/* calculate the Eij coefficients */
					CalcEijt(Eijx);
					CalcEijt(Eijy);
					CalcEijt(Eijz);

					Ex0_1 = Ex0 * exp(-ai * aj * ABx * ABx / aij);
					Ey0_1 = Ey0 * exp(-ai * aj * ABy * ABy / aij);
					Ez0_1 = Ez0 * exp(-ai * aj * ABz * ABz / aij);
					Exyz = Ex0_1 * Ey0_1 * Ez0_1;

					kPP = kPx * Px + kPy * Py + kPz * Pz;
					fac = pow(M_PI / aij, 1.5) * exp(-0.25 * kP2 / aij) * exp(I * kPP);

					/* calculate R_tuv integrals - overlap/kinetic */
					ROvrl.zero();
					ROvrl.load(kPx, kPy, kPz, Px, Py, Pz, aij);
					CalcROvrl(ROvrl, fac);

					/* radial CAP integrals */
					cap_rad.zero();
					cap_tmp.zero();
					if (keys.capi) {
						CAPRadFac(lij, lij + 4, sqrt(kP2) * keys.cap_on,
						          aij * keys.cap_on * keys.cap_on, cap_rad.v);

						///for(int l=0; l<=lij; l++) for(int n=l+2; n<=lij+4; n++)
						///    cout << l << " " << n << " " << cap_rad.v[l][n] << endl;

						usint mij;
						double pre_fac;
						pre_fac = 4.0 * M_PI * pow(keys.cap_on, lij + 5);

						for (mij = 0; mij < cap_tmp.num; mij++) {
							II = 1.0;
							for (usint L = 0; L <= lij; L++) {
								cap_tmp[mij] += II * cap_ang.v[L][mij] *
								                (cap_rad.v[L][lij + 4] - 2.0 * cap_rad.v[L][lij + 3] + cap_rad.v[L][lij + 2]);
								II = II * I;
							};
							cap_tmp[mij] = cap_tmp[mij] * pre_fac;
						};
					};

					ia = 0;
					ja = 0;
					/* loop over Cartesian components - shell I */
					for (mi = 0; mi < shgi; mi++) {
						ka = li - ia - ja;
						ib = 0;
						jb = 0;
						/* loop over Cartesian components - shell J */
						for (mj = 0; mj < shgj; mj++) {
							kb = lj - ib - jb;

							/* filling the hash table */
							for (int hi = 0; hi < dim_hash; hi++) {
								ia1 = ia + hash[hi][0];
								ja1 = ja + hash[hi][1];
								ka1 = ka + hash[hi][2];
								if (ia1 < 0 || ja1 < 0 || ka1 < 0)
									continue;

								for (int hj = 0; hj < dim_hash; hj++) {
									ib1 = ib + hash[hj][0];
									jb1 = jb + hash[hj][1];
									kb1 = kb + hash[hj][2];
									if (ib1 < 0 || jb1 < 0 || kb1 < 0)
										continue;

									tmp3s = 0.0;
									for (int t = 0; t <= ia1 + ib1; t++) {
										tmp2s = 0.0;
										for (int u = 0; u <= ja1 + jb1; u++) {
											tmp1s = 0.0;
											for (int v = 0; v <= ka1 + kb1; v++) {
												tmp1s += ROvrl.rtuv[t][u][v] * Eijz.v[ka1][kb1][v];
											};
											tmp2s += tmp1s * Eijy.v[ja1][jb1][u];
										};
										tmp3s += tmp2s * Eijx.v[ia1][ib1][t];
									};

									hash_fill[hi][hj] = tmp3s;
								};
							};

							/* overlap integrals */
							shl_crt_s[mi * shgj + mj] += hash_fill[6][6] * dij * Exyz;

							/* kinetic energy integrals */
							tmp1s = 0.0;
							tmp2s = 0.0;
							tmp3s = 0.0;

							tmp1s = (hash_fill[0][0] * ib - 2.0 * aj * hash_fill[0][3] + I * basis[j].kx * hash_fill[0][6]) * ia - (hash_fill[3][0] * ib - 2.0 * aj * hash_fill[3][3] + I * basis[j].kx * hash_fill[3][6]) * 2.0 * ai - (hash_fill[6][0] * ib - 2.0 * aj * hash_fill[6][3] + I * basis[j].kx * hash_fill[6][6]) * I * basis[i].kx;

							tmp2s = (hash_fill[1][1] * jb - 2.0 * aj * hash_fill[1][4] + I * basis[j].ky * hash_fill[1][6]) * ja - (hash_fill[4][1] * jb - 2.0 * aj * hash_fill[4][4] + I * basis[j].ky * hash_fill[4][6]) * 2.0 * ai - (hash_fill[6][1] * jb - 2.0 * aj * hash_fill[6][4] + I * basis[j].ky * hash_fill[6][6]) * I * basis[i].ky;

							tmp3s = (hash_fill[2][2] * kb - 2.0 * aj * hash_fill[2][5] + I * basis[j].kz * hash_fill[2][6]) * ka - (hash_fill[5][2] * kb - 2.0 * aj * hash_fill[5][5] + I * basis[j].kz * hash_fill[5][6]) * 2.0 * ai - (hash_fill[6][2] * kb - 2.0 * aj * hash_fill[6][5] + I * basis[j].kz * hash_fill[6][6]) * I * basis[i].kz;

							shl_crt_t[mi * shgj + mj] += (tmp1s + tmp2s + tmp3s) * dij * Exyz / 2.0;

							shl_crt_tx[mi * shgj + mj] += tmp1s * dij * Exyz / 2.0;
							shl_crt_ty[mi * shgj + mj] += tmp2s * dij * Exyz / 2.0;
							shl_crt_tz[mi * shgj + mj] += tmp3s * dij * Exyz / 2.0;

							/* dipole moment integrals */
							if (keys.dip) {
								BPx = basis[j].Ax - std::get<0>(keys.points[0]);
								BPy = basis[j].Ay - std::get<1>(keys.points[0]);
								BPz = basis[j].Az - std::get<2>(keys.points[0]);

								tmp1s = hash_fill[6][3] + BPx * hash_fill[6][6];
								tmp2s = hash_fill[6][4] + BPy * hash_fill[6][6];
								tmp3s = hash_fill[6][5] + BPz * hash_fill[6][6];

								shl_crt_dx[mi * shgj + mj] += tmp1s * dij * Exyz;
								shl_crt_dy[mi * shgj + mj] += tmp2s * dij * Exyz;
								shl_crt_dz[mi * shgj + mj] += tmp3s * dij * Exyz;
							};

							/* quadruple moment integrals */
							if (keys.quad) {
								APx = basis[i].Ax - std::get<0>(keys.points[0]);
								APy = basis[i].Ay - std::get<1>(keys.points[0]);
								APz = basis[i].Az - std::get<2>(keys.points[0]);

								BPx = basis[j].Ax - std::get<0>(keys.points[0]);
								BPy = basis[j].Ay - std::get<1>(keys.points[0]);
								BPz = basis[j].Az - std::get<2>(keys.points[0]);

								/* x-x */
								tmp1s = hash_fill[3][3] + hash_fill[3][6] * BPx + hash_fill[6][3] * APx + hash_fill[6][6] * APx * BPx;

								shl_crt_xx[mi * shgj + mj] += tmp1s * dij * Exyz;

								/* x-y */
								tmp1s = hash_fill[3][4] + hash_fill[3][6] * BPy + hash_fill[6][4] * APx + hash_fill[6][6] * APx * BPy;

								shl_crt_xy[mi * shgj + mj] += tmp1s * dij * Exyz;

								/* x-z */
								tmp1s = hash_fill[3][5] + hash_fill[3][6] * BPz + hash_fill[6][5] * APx + hash_fill[6][6] * APx * BPz;

								shl_crt_xz[mi * shgj + mj] += tmp1s * dij * Exyz;

								/* y-y */
								tmp1s = hash_fill[4][4] + hash_fill[4][6] * BPy + hash_fill[6][4] * APy + hash_fill[6][6] * APy * BPy;

								shl_crt_yy[mi * shgj + mj] += tmp1s * dij * Exyz;

								/* y-z */
								tmp1s = hash_fill[4][5] + hash_fill[4][6] * BPz + hash_fill[6][5] * APy + hash_fill[6][6] * APy * BPz;

								shl_crt_yz[mi * shgj + mj] += tmp1s * dij * Exyz;

								/* z-z */
								tmp1s = hash_fill[5][5] + hash_fill[5][6] * BPz + hash_fill[6][5] * APz + hash_fill[6][6] * APz * BPz;

								shl_crt_zz[mi * shgj + mj] += tmp1s * dij * Exyz;
							};

							/* dipole velocity integrals */
							if (keys.velo) {
								tmp1s = hash_fill[6][0] * ib - 2.0 * aj * hash_fill[6][3] + I * basis[j].kx * hash_fill[6][6];
								tmp2s = hash_fill[6][1] * jb - 2.0 * aj * hash_fill[6][4] + I * basis[j].ky * hash_fill[6][6];
								tmp3s = hash_fill[6][2] * kb - 2.0 * aj * hash_fill[6][5] + I * basis[j].kz * hash_fill[6][6];

								shl_crt_gx[mi * shgj + mj] += tmp1s * dij * Exyz;
								shl_crt_gy[mi * shgj + mj] += tmp2s * dij * Exyz;
								shl_crt_gz[mi * shgj + mj] += tmp3s * dij * Exyz;
							};

							/* angular momentum operator, L_z */

							/* complex absorbing potentials integrals */
							if (keys.capi) {
								usint pos;
								pos = FindPos(ia + ib, ja + jb, ka + kb);
								shl_crt_ca[mi * shgj + mj] += dij * keys.cap_mul * cap_tmp[pos];
							};

							jb++;
							if (jb > lj - ib) {
								jb = 0;
								ib++;
							};
						};
						ja++;
						if (ja > li - ia) {
							ja = 0;
							ia++;
						};
					};

					/* loop over all nuclei */
					for (usint nc = 0; nc < nucli.size(); nc++) {
						/* calculate R_tuv integrals - nuclear attr. */
						RNucA.zero();
						RNucA.load(kPx, kPy, kPz, Px, Py, Pz,
						           nucli[nc].Cx, nucli[nc].Cy,
						           nucli[nc].Cz, aij);
						CalcRNucA(RNucA);

						ia = 0;
						ja = 0;
						/* loop over Cartesian components - shell I */
						for (mi = 0; mi < shgi; mi++) {
							ka = li - ia - ja;
							ib = 0;
							jb = 0;
							/* loop over Cartesian components - shell J */
							for (mj = 0; mj < shgj; mj++) {
								kb = lj - ib - jb;

								tmp3s = 0.0;
								for (int t = 0; t <= ia + ib; t++) {
									tmp2s = 0.0;
									for (int u = 0; u <= ja + jb; u++) {
										tmp1s = 0.0;
										for (int v = 0; v <= ka + kb; v++) {
											tmp1s += RNucA.rtuv[t][u][v] * Eijz.v[ka][kb][v];
										};
										tmp2s += tmp1s * Eijy.v[ja][jb][u];
									};
									tmp3s += tmp2s * Eijx.v[ia][ib][t];
								};

								shl_crt_v[mi * shgj + mj] -= tmp3s * dij * Exyz * nucli[nc].chrg;

								jb++;
								if (jb > lj - ib) {
									jb = 0;
									ib++;
								};
							};
							ja++;
							if (ja > li - ia) {
								ja = 0;
								ia++;
							};
						};

					}; /* end loop over nuclei */
				};
			}; /* end loop over contractions */

			/* rescale to fix the angular normalisation */
			ia = 0;
			ja = 0;
			for (mi = 0; mi < shgi; mi++) {
				ka = li - ia - ja;
				ib = 0;
				jb = 0;
				for (mj = 0; mj < shgj; mj++) {
					kb = lj - ib - jb;
					shl_crt_s[mi * shgj + mj] = shl_crt_s[mi * shgj + mj] / xyz_norm[li][mi] / xyz_norm[lj][mj];
					shl_crt_v[mi * shgj + mj] = shl_crt_v[mi * shgj + mj] / xyz_norm[li][mi] / xyz_norm[lj][mj];
					shl_crt_t[mi * shgj + mj] = shl_crt_t[mi * shgj + mj] / xyz_norm[li][mi] / xyz_norm[lj][mj];

					if (keys.dip) {
						shl_crt_dx[mi * shgj + mj] = shl_crt_dx[mi * shgj + mj] / xyz_norm[li][mi] / xyz_norm[lj][mj];
						shl_crt_dy[mi * shgj + mj] = shl_crt_dy[mi * shgj + mj] / xyz_norm[li][mi] / xyz_norm[lj][mj];
						shl_crt_dz[mi * shgj + mj] = shl_crt_dz[mi * shgj + mj] / xyz_norm[li][mi] / xyz_norm[lj][mj];
					};

					if (keys.quad) {
						shl_crt_xx[mi * shgj + mj] = shl_crt_xx[mi * shgj + mj] / xyz_norm[li][mi] / xyz_norm[lj][mj];
						shl_crt_yy[mi * shgj + mj] = shl_crt_yy[mi * shgj + mj] / xyz_norm[li][mi] / xyz_norm[lj][mj];
						shl_crt_zz[mi * shgj + mj] = shl_crt_zz[mi * shgj + mj] / xyz_norm[li][mi] / xyz_norm[lj][mj];
						shl_crt_xy[mi * shgj + mj] = shl_crt_xy[mi * shgj + mj] / xyz_norm[li][mi] / xyz_norm[lj][mj];
						shl_crt_xz[mi * shgj + mj] = shl_crt_xz[mi * shgj + mj] / xyz_norm[li][mi] / xyz_norm[lj][mj];
						shl_crt_yz[mi * shgj + mj] = shl_crt_yz[mi * shgj + mj] / xyz_norm[li][mi] / xyz_norm[lj][mj];
					};

					if (keys.velo) {
						shl_crt_gx[mi * shgj + mj] = shl_crt_gx[mi * shgj + mj] / xyz_norm[li][mi] / xyz_norm[lj][mj];
						shl_crt_gy[mi * shgj + mj] = shl_crt_gy[mi * shgj + mj] / xyz_norm[li][mi] / xyz_norm[lj][mj];
						shl_crt_gz[mi * shgj + mj] = shl_crt_gz[mi * shgj + mj] / xyz_norm[li][mi] / xyz_norm[lj][mj];
					};

					if (keys.kinc) {
						shl_crt_tx[mi * shgj + mj] = shl_crt_tx[mi * shgj + mj] / xyz_norm[li][mi] / xyz_norm[lj][mj];
						shl_crt_ty[mi * shgj + mj] = shl_crt_ty[mi * shgj + mj] / xyz_norm[li][mi] / xyz_norm[lj][mj];
						shl_crt_tz[mi * shgj + mj] = shl_crt_tz[mi * shgj + mj] / xyz_norm[li][mi] / xyz_norm[lj][mj];
					};

					if (keys.capi)
						shl_crt_ca[mi * shgj + mj] = shl_crt_ca[mi * shgj + mj] / xyz_norm[li][mi] / xyz_norm[lj][mj];

					jb++;
					if (jb > lj - ib) {
						jb = 0;
						ib++;
					};
				};
				ja++;
				if (ja > li - ia) {
					ja = 0;
					ia++;
				};
			};

			/* transform to the spherical representation */
			TransToSpher(shl_crt_s.v, shl_sph_s.v, li, lj);
			TransToSpher(shl_crt_t.v, shl_sph_t.v, li, lj);
			TransToSpher(shl_crt_v.v, shl_sph_v.v, li, lj);

			if (keys.dip) {
				TransToSpher(shl_crt_dx.v, shl_sph_dx.v, li, lj);
				TransToSpher(shl_crt_dy.v, shl_sph_dy.v, li, lj);
				TransToSpher(shl_crt_dz.v, shl_sph_dz.v, li, lj);
			};

			if (keys.quad) {
				TransToSpher(shl_crt_xx.v, shl_sph_xx.v, li, lj);
				TransToSpher(shl_crt_yy.v, shl_sph_yy.v, li, lj);
				TransToSpher(shl_crt_zz.v, shl_sph_zz.v, li, lj);
				TransToSpher(shl_crt_xy.v, shl_sph_xy.v, li, lj);
				TransToSpher(shl_crt_xz.v, shl_sph_xz.v, li, lj);
				TransToSpher(shl_crt_yz.v, shl_sph_yz.v, li, lj);
			};

			if (keys.velo) {
				TransToSpher(shl_crt_gx.v, shl_sph_gx.v, li, lj);
				TransToSpher(shl_crt_gy.v, shl_sph_gy.v, li, lj);
				TransToSpher(shl_crt_gz.v, shl_sph_gz.v, li, lj);
			};

			if (keys.kinc) {
				TransToSpher(shl_crt_tx.v, shl_sph_tx.v, li, lj);
				TransToSpher(shl_crt_ty.v, shl_sph_ty.v, li, lj);
				TransToSpher(shl_crt_tz.v, shl_sph_tz.v, li, lj);
			};

			if (keys.capi)
				TransToSpher(shl_crt_ca.v, shl_sph_ca.v, li, lj);

			/* switch to the Gamess holy order */
			if (gamess_order) {
				shl_crt_dum.zero();
				GamessHolyOrder(shl_crt_s.v, shl_crt_dum.v, li, lj);
				memcpy(shl_crt_s.v, shl_crt_dum.v, shl_siz_crt2 * sizeof(cdouble));

				shl_crt_dum.zero();
				GamessHolyOrder(shl_crt_t.v, shl_crt_dum.v, li, lj);
				memcpy(shl_crt_t.v, shl_crt_dum.v, shl_siz_crt2 * sizeof(cdouble));

				shl_crt_dum.zero();
				GamessHolyOrder(shl_crt_v.v, shl_crt_dum.v, li, lj);
				memcpy(shl_crt_v.v, shl_crt_dum.v, shl_siz_crt2 * sizeof(cdouble));

				if (keys.dip) {
					shl_crt_dum.zero();
					GamessHolyOrder(shl_crt_dx.v, shl_crt_dum.v, li, lj);
					memcpy(shl_crt_dx.v, shl_crt_dum.v, shl_siz_crt2 * sizeof(cdouble));

					shl_crt_dum.zero();
					GamessHolyOrder(shl_crt_dy.v, shl_crt_dum.v, li, lj);
					memcpy(shl_crt_dy.v, shl_crt_dum.v, shl_siz_crt2 * sizeof(cdouble));

					shl_crt_dum.zero();
					GamessHolyOrder(shl_crt_dz.v, shl_crt_dum.v, li, lj);
					memcpy(shl_crt_dz.v, shl_crt_dum.v, shl_siz_crt2 * sizeof(cdouble));
				};

				if (keys.quad) {
					shl_crt_dum.zero();
					GamessHolyOrder(shl_crt_xx.v, shl_crt_dum.v, li, lj);
					memcpy(shl_crt_xx.v, shl_crt_dum.v, shl_siz_crt2 * sizeof(cdouble));

					shl_crt_dum.zero();
					GamessHolyOrder(shl_crt_yy.v, shl_crt_dum.v, li, lj);
					memcpy(shl_crt_yy.v, shl_crt_dum.v, shl_siz_crt2 * sizeof(cdouble));

					shl_crt_dum.zero();
					GamessHolyOrder(shl_crt_zz.v, shl_crt_dum.v, li, lj);
					memcpy(shl_crt_zz.v, shl_crt_dum.v, shl_siz_crt2 * sizeof(cdouble));

					shl_crt_dum.zero();
					GamessHolyOrder(shl_crt_xy.v, shl_crt_dum.v, li, lj);
					memcpy(shl_crt_xy.v, shl_crt_dum.v, shl_siz_crt2 * sizeof(cdouble));

					shl_crt_dum.zero();
					GamessHolyOrder(shl_crt_xz.v, shl_crt_dum.v, li, lj);
					memcpy(shl_crt_xz.v, shl_crt_dum.v, shl_siz_crt2 * sizeof(cdouble));

					shl_crt_dum.zero();
					GamessHolyOrder(shl_crt_yz.v, shl_crt_dum.v, li, lj);
					memcpy(shl_crt_yz.v, shl_crt_dum.v, shl_siz_crt2 * sizeof(cdouble));
				};

				if (keys.velo) {
					shl_crt_dum.zero();
					GamessHolyOrder(shl_crt_gx.v, shl_crt_dum.v, li, lj);
					memcpy(shl_crt_gx.v, shl_crt_dum.v, shl_siz_crt2 * sizeof(cdouble));

					shl_crt_dum.zero();
					GamessHolyOrder(shl_crt_gy.v, shl_crt_dum.v, li, lj);
					memcpy(shl_crt_gy.v, shl_crt_dum.v, shl_siz_crt2 * sizeof(cdouble));

					shl_crt_dum.zero();
					GamessHolyOrder(shl_crt_gz.v, shl_crt_dum.v, li, lj);
					memcpy(shl_crt_gz.v, shl_crt_dum.v, shl_siz_crt2 * sizeof(cdouble));
				};

				if (keys.kinc) {
					shl_crt_dum.zero();
					GamessHolyOrder(shl_crt_tx.v, shl_crt_dum.v, li, lj);
					memcpy(shl_crt_tx.v, shl_crt_dum.v, shl_siz_crt2 * sizeof(cdouble));

					shl_crt_dum.zero();
					GamessHolyOrder(shl_crt_ty.v, shl_crt_dum.v, li, lj);
					memcpy(shl_crt_ty.v, shl_crt_dum.v, shl_siz_crt2 * sizeof(cdouble));

					shl_crt_dum.zero();
					GamessHolyOrder(shl_crt_tz.v, shl_crt_dum.v, li, lj);
					memcpy(shl_crt_tz.v, shl_crt_dum.v, shl_siz_crt2 * sizeof(cdouble));
				};

				if (keys.capi) {
					shl_crt_dum.zero();
					GamessHolyOrder(shl_crt_ca.v, shl_crt_dum.v, li, lj);
					memcpy(shl_crt_ca.v, shl_crt_dum.v, shl_siz_crt2 * sizeof(cdouble));
				};
			};

			/* place in the table - Cartesian rep. */
			for (mi = 0; mi < shgi; mi++) {
				indi = posimx_crt + mi;
				mjmx = (i == j) ? mi + 1 : shgj;
				for (mj = 0; mj < mjmx; mj++) {
					indj = posjmx_crt + mj;

					ovrl_crt[indi * num_crt + indj] = shl_crt_s[mi * shgj + mj];
					ovrl_crt[indj * num_crt + indi] = conj(shl_crt_s[mi * shgj + mj]);

					nuc_crt[indi * num_crt + indj] = shl_crt_v[mi * shgj + mj];
					nuc_crt[indj * num_crt + indi] = conj(shl_crt_v[mi * shgj + mj]);

					kin_crt[indi * num_crt + indj] = shl_crt_t[mi * shgj + mj];
					kin_crt[indj * num_crt + indi] = conj(shl_crt_t[mi * shgj + mj]);

					if (keys.dip) {
						dipx_crt[indi * num_crt + indj] = shl_crt_dx[mi * shgj + mj];
						dipx_crt[indj * num_crt + indi] = conj(shl_crt_dx[mi * shgj + mj]);

						dipy_crt[indi * num_crt + indj] = shl_crt_dy[mi * shgj + mj];
						dipy_crt[indj * num_crt + indi] = conj(shl_crt_dy[mi * shgj + mj]);

						dipz_crt[indi * num_crt + indj] = shl_crt_dz[mi * shgj + mj];
						dipz_crt[indj * num_crt + indi] = conj(shl_crt_dz[mi * shgj + mj]);
					};

					if (keys.quad) {
						qdxx_crt[indi * num_crt + indj] = shl_crt_xx[mi * shgj + mj];
						qdxx_crt[indj * num_crt + indi] = conj(shl_crt_xx[mi * shgj + mj]);

						qdyy_crt[indi * num_crt + indj] = shl_crt_yy[mi * shgj + mj];
						qdyy_crt[indj * num_crt + indi] = conj(shl_crt_yy[mi * shgj + mj]);

						qdzz_crt[indi * num_crt + indj] = shl_crt_zz[mi * shgj + mj];
						qdzz_crt[indj * num_crt + indi] = conj(shl_crt_zz[mi * shgj + mj]);

						qdxy_crt[indi * num_crt + indj] = shl_crt_xy[mi * shgj + mj];
						qdxy_crt[indj * num_crt + indi] = conj(shl_crt_xy[mi * shgj + mj]);

						qdxz_crt[indi * num_crt + indj] = shl_crt_xz[mi * shgj + mj];
						qdxz_crt[indj * num_crt + indi] = conj(shl_crt_xz[mi * shgj + mj]);

						qdyz_crt[indi * num_crt + indj] = shl_crt_yz[mi * shgj + mj];
						qdyz_crt[indj * num_crt + indi] = conj(shl_crt_yz[mi * shgj + mj]);
					};

					if (keys.velo) {
						grdx_crt[indi * num_crt + indj] = shl_crt_gx[mi * shgj + mj];
						grdx_crt[indj * num_crt + indi] = -conj(shl_crt_gx[mi * shgj + mj]);

						grdy_crt[indi * num_crt + indj] = shl_crt_gy[mi * shgj + mj];
						grdy_crt[indj * num_crt + indi] = -conj(shl_crt_gy[mi * shgj + mj]);

						grdz_crt[indi * num_crt + indj] = shl_crt_gz[mi * shgj + mj];
						grdz_crt[indj * num_crt + indi] = -conj(shl_crt_gz[mi * shgj + mj]);
					};

					if (keys.kinc) {
						tx_crt[indi * num_crt + indj] = shl_crt_tx[mi * shgj + mj];
						tx_crt[indj * num_crt + indi] = conj(shl_crt_tx[mi * shgj + mj]);

						ty_crt[indi * num_crt + indj] = shl_crt_ty[mi * shgj + mj];
						ty_crt[indj * num_crt + indi] = conj(shl_crt_ty[mi * shgj + mj]);

						tz_crt[indi * num_crt + indj] = shl_crt_tz[mi * shgj + mj];
						tz_crt[indj * num_crt + indi] = conj(shl_crt_tz[mi * shgj + mj]);
					};

					if (keys.capi) {
						capi_crt[indi * num_crt + indj] = shl_crt_ca[mi * shgj + mj];
						capi_crt[indj * num_crt + indi] = conj(shl_crt_ca[mi * shgj + mj]);
					};
				};
			};

			/* place in the table - spherical rep. */
			for (mi = 0; mi < shi; mi++) {
				indi = posimx_sph + mi;
				mjmx = (i == j) ? mi + 1 : shj;
				for (mj = 0; mj < mjmx; mj++) {
					indj = posjmx_sph + mj;

					ovrl_sph[indi * num_sph + indj] = shl_sph_s[mi * shj + mj];
					ovrl_sph[indj * num_sph + indi] = conj(shl_sph_s[mi * shj + mj]);

					nuc_sph[indi * num_sph + indj] = shl_sph_v[mi * shj + mj];
					nuc_sph[indj * num_sph + indi] = conj(shl_sph_v[mi * shj + mj]);

					kin_sph[indi * num_sph + indj] = shl_sph_t[mi * shj + mj];
					kin_sph[indj * num_sph + indi] = conj(shl_sph_t[mi * shj + mj]);

					if (keys.dip) {
						dipx_sph[indi * num_sph + indj] = shl_sph_dx[mi * shj + mj];
						dipx_sph[indj * num_sph + indi] = conj(shl_sph_dx[mi * shj + mj]);

						dipy_sph[indi * num_sph + indj] = shl_sph_dy[mi * shj + mj];
						dipy_sph[indj * num_sph + indi] = conj(shl_sph_dy[mi * shj + mj]);

						dipz_sph[indi * num_sph + indj] = shl_sph_dz[mi * shj + mj];
						dipz_sph[indj * num_sph + indi] = conj(shl_sph_dz[mi * shj + mj]);
					};

					if (keys.quad) {
						qdxx_sph[indi * num_sph + indj] = shl_sph_xx[mi * shj + mj];
						qdxx_sph[indj * num_sph + indi] = conj(shl_sph_xx[mi * shj + mj]);

						qdyy_sph[indi * num_sph + indj] = shl_sph_yy[mi * shj + mj];
						qdyy_sph[indj * num_sph + indi] = conj(shl_sph_yy[mi * shj + mj]);

						qdzz_sph[indi * num_sph + indj] = shl_sph_zz[mi * shj + mj];
						qdzz_sph[indj * num_sph + indi] = conj(shl_sph_zz[mi * shj + mj]);

						qdxy_sph[indi * num_sph + indj] = shl_sph_xy[mi * shj + mj];
						qdxy_sph[indj * num_sph + indi] = conj(shl_sph_xy[mi * shj + mj]);

						qdxz_sph[indi * num_sph + indj] = shl_sph_xz[mi * shj + mj];
						qdxz_sph[indj * num_sph + indi] = conj(shl_sph_xz[mi * shj + mj]);

						qdyz_sph[indi * num_sph + indj] = shl_sph_yz[mi * shj + mj];
						qdyz_sph[indj * num_sph + indi] = conj(shl_sph_yz[mi * shj + mj]);
					};

					if (keys.velo) {
						grdx_sph[indi * num_sph + indj] = shl_sph_gx[mi * shj + mj];
						grdx_sph[indj * num_sph + indi] = -conj(shl_sph_gx[mi * shj + mj]);

						grdy_sph[indi * num_sph + indj] = shl_sph_gy[mi * shj + mj];
						grdy_sph[indj * num_sph + indi] = -conj(shl_sph_gy[mi * shj + mj]);

						grdz_sph[indi * num_sph + indj] = shl_sph_gz[mi * shj + mj];
						grdz_sph[indj * num_sph + indi] = -conj(shl_sph_gz[mi * shj + mj]);
					};

					if (keys.kinc) {
						tx_sph[indi * num_sph + indj] = shl_sph_tx[mi * shj + mj];
						tx_sph[indj * num_sph + indi] = conj(shl_sph_tx[mi * shj + mj]);

						ty_sph[indi * num_sph + indj] = shl_sph_ty[mi * shj + mj];
						ty_sph[indj * num_sph + indi] = conj(shl_sph_ty[mi * shj + mj]);

						tz_sph[indi * num_sph + indj] = shl_sph_tz[mi * shj + mj];
						tz_sph[indj * num_sph + indi] = conj(shl_sph_tz[mi * shj + mj]);
					};

					if (keys.capi) {
						capi_sph[indi * num_sph + indj] = shl_sph_ca[mi * shj + mj];
						capi_sph[indj * num_sph + indi] = conj(shl_sph_ca[mi * shj + mj]);
					};
				};
			};

			posjmx_crt += shgj;
			posjmx_sph += shj;
		};
		posimx_crt += shgi;
		posimx_sph += shi;
	};

	for (int i = 0; i < num_crt2; i++)
		bare_crt[i] = kin_crt[i] + nuc_crt[i];
	for (int i = 0; i < num_sph2; i++)
		bare_sph[i] = kin_sph[i] + nuc_sph[i];

	///PrintCMatrix( nuc_crt.v , num_crt );
	///PrintCMatrix( ovrl_crt.v, num_crt );
	///PrintCMatrix( kin_crt.v , num_crt );
	///PrintCMatrix( bare_crt.v, num_crt );
	///PrintCMatrix( dipz_crt.v, num_crt );
	///PrintCMatrix( ovrl_crt.v, num_crt );
	///PrintCMatrix( capi_crt.v, num_crt );

	//PrintCMatrix( ovrl_crt.v, num_crt );

	/* write the one-electron integrals to the disk - cartesian */
	std::ofstream ofs(keys.file1E, std::ios::out | std::ios::binary);

	/*check the file1E path (MS)*/
	if (!ofs.is_open()) {
		cout << "Cannot open file1E!\n";
		return EXIT_FAILURE;
	}

	WriteDown(ovrl_crt.v, num_crt2, ofs);
	WriteDown(kin_crt.v, num_crt2, ofs);
	WriteDown(nuc_crt.v, num_crt2, ofs);
	WriteDown(bare_crt.v, num_crt2, ofs);

	WriteDown(dipx_crt.v, num_crt2, ofs);
	WriteDown(dipy_crt.v, num_crt2, ofs);
	WriteDown(dipz_crt.v, num_crt2, ofs);

	WriteDown(qdxx_crt.v, num_crt2, ofs);
	WriteDown(qdyy_crt.v, num_crt2, ofs);
	WriteDown(qdzz_crt.v, num_crt2, ofs);
	WriteDown(qdxy_crt.v, num_crt2, ofs);
	WriteDown(qdxz_crt.v, num_crt2, ofs);
	WriteDown(qdyz_crt.v, num_crt2, ofs);

	WriteDown(grdx_crt.v, num_crt2, ofs);
	WriteDown(grdy_crt.v, num_crt2, ofs);
	WriteDown(grdz_crt.v, num_crt2, ofs);

	WriteDown(tx_crt.v, num_crt2, ofs);
	WriteDown(ty_crt.v, num_crt2, ofs);
	WriteDown(tz_crt.v, num_crt2, ofs);

	WriteDown(capi_crt.v, num_crt2, ofs);

	ofs.close();

	/* write the one-electron integrals to the disk - spherical */
	/*
=======
int main(int argc, char* argv[]) {
    /* put a battery in the clock */
    double t_start, t_end;
    t_start = omp_get_wtime();
    
    /* code initialisation */
    Keywords keys;
    Init();
    string inpname;
    if( argc > 1 ) { inpname = argv[1]; } else
        { cout << " Name of the input file not specified! Emergency halt." << endl; exit( EXIT_FAILURE ); };
    keys.name_me( inpname );
    
    /* read the $BASIS keyword */
    vector <GaussC> basis; 
    vector <Nuclei> nucli;
    
    GaussC gau_tmp;
    Nuclei nuc_tmp;
    
    int num_crt,num_sph;
    num_crt = 0;
    num_sph = 0;
    
    ifstream ifile;
    string line,name;
    stringstream ss;
    ifile.open ( inpname );
    
    do { getline(ifile, line); } while( line != "$BASIS" );
    getline(ifile, line);
    do {
        ss.clear();
        ss << line;
        line.clear();
        ss >> nuc_tmp.name >> nuc_tmp.chrg >> nuc_tmp.Cx >> nuc_tmp.Cy >> nuc_tmp.Cz;
        
        nucli.push_back( std::move( nuc_tmp ) );
        getline( ifile, line );
        string dummy, moment;
        int end;
    
        double ex,d_re,d_im;
        do {
            ss.clear();
            ss << line ;
            line.clear();
            ss >> moment >> end;
            
            if(moment == "S")      gau_tmp.lA = 0;
            else if(moment == "P") gau_tmp.lA = 1;
            else if(moment == "D") gau_tmp.lA = 2;
            else if(moment == "F") gau_tmp.lA = 3;
            else if(moment == "G") gau_tmp.lA = 4;
            else if(moment == "H") gau_tmp.lA = 5;
            else if(moment == "I") gau_tmp.lA = 6;
            else if(moment == "K") gau_tmp.lA = 7;
            else if(moment == "L") gau_tmp.lA = 8;
            else {
                cout << " Unrecognised angular momentum of a function!" << endl;
                cout << " Emergency halt." << endl;
                exit( EXIT_FAILURE );
            };
            
            num_crt += crt_siz[ gau_tmp.lA ];
            num_sph += sph_siz[ gau_tmp.lA ];
            
            gau_tmp.clen = end;
            for(int i=0; i<end; i++) {
                getline( ifile, line );
                ss.clear();
                ss << line;
                line.clear();
                ss >> dummy >> ex >> d_re >> d_im >> gau_tmp.kx >> gau_tmp.ky >> gau_tmp.kz;
                gau_tmp.alphaA.push_back( ex );
                gau_tmp.dA_re.push_back( d_re );
                gau_tmp.dA_im.push_back( d_im );
            };
            
            RenormContr( gau_tmp, keys.norm1E );
            gau_tmp.Ax = nuc_tmp.Cx;
            gau_tmp.Ay = nuc_tmp.Cy;
            gau_tmp.Az = nuc_tmp.Cz;
    
            basis.push_back( std::move( gau_tmp ) );
            getline( ifile, line );
        } while ( !line.empty() );
        getline(ifile,line);
    } while( line != "$END" );
    ifile.close();
    
    /* basis set data - mostly optional printing */
    cout << endl;
    cout << " Number of basis set shells      = " << basis.size() << endl;
    cout << " Number of functions (spherical) = " << num_sph << endl;
    cout << " Number of functions (cartesian) = " << num_crt << endl;
    cout << " Contracted functions are now normalised to the unity." << endl;
    cout << endl;
    
    /* read the remaining keywords */
    ifile.open ( inpname );
    while( getline(ifile, line) ) if( line == "$INTS" ) break;
    getline(ifile, line);
    int num_ints = atoi( line.c_str() );
    
    for(int i=0; i<num_ints; i++) {
        getline(ifile, line);
        if( line == "STVH" ) keys.stvh = true;
        else if ( line == "DIPOLE"   ) keys.dip    = true;
        else if ( line == "QUADRUP"  ) keys.quad   = true;
        else if ( line == "VELOCITY" ) keys.velo   = true;
        else if ( line == "ERI"      ) keys.eri    = true;
        else if ( line == "PROJ"     ) keys.proj   = true;
        else if ( line == "KINXYZ"   ) keys.kinc   = true;
        else if ( line == "CAPINT"   ) keys.capi   = true;
        else cout << " Unrecognised integral type (" << i + 1 << "). Omitted." << endl;
    };
    
    cout << " Requested one-electron integrals:" << endl;
    if( keys.stvh ) cout << " - ordinary one-electron STVH integrals" << endl;
    if( keys.dip  ) cout << " - dipole moment integrals" << endl;
    if( keys.quad ) cout << " - quadruple moment integrals" << endl;
    if( keys.velo ) cout << " - dipole velocity integrals" << endl;
    if( keys.proj ) cout << " - projection integrals" << endl;
    if( keys.kinc ) cout << " - kinetic energy components int." << endl;
    if( keys.capi ) cout << " - complex absorbing potentials int." << endl;
    cout << endl;
    ifile.close();
    
    int num_points(0);
    line.clear();
    ifile.open ( inpname );
    while( getline(ifile, line) ) if( line == "$POINTS" ) {
        getline(ifile, line);
        num_points = atoi( line.c_str() );
        break;
    };
    
    cout << " Points for properties calculations (" << num_points << "): " << endl;
    for(int i=0; i<num_points; i++) {
        getline(ifile, line);
        ss.clear();
        ss << line ;
        double xp,yp,zp;
        ss >> xp >> yp >> zp;
        keys.points.push_back( std::make_tuple( xp, yp, zp ) );
        cout << " " << std::get<0>( keys.points[i] ) << " " << std::get<1>( keys.points[i] ) << " " << std::get<2>( keys.points[i] ) << endl;
    };
    ifile.close();
    
    line.clear();
    ifile.open ( inpname );
    while( getline(ifile, line) ) if( line == "$RGRID" ) { 
        getline(ifile, line); keys.n_rad     = atoi( line.c_str() );
        getline(ifile, line); keys.rad_str   = atof( line.c_str() );
        getline(ifile, line); keys.rad_stp   = atof( line.c_str() );
        getline(ifile, line); keys.lproj_max = atoi( line.c_str() );
        break;
    };
    
    cout << endl;
    if( keys.n_rad     <= 0   || 
        keys.rad_stp   <= 0.0 || 
        keys.rad_str   <  0.0  )
        keys.proj = false; /// fail-safe
    
    if( keys.proj ) {
        cout << " Radial grid parameters:" << endl;
        cout << "  Number of radial points     = " << keys.n_rad << endl;
        cout << "  Maximal projection momentum = " << keys.lproj_max << endl;
        cout << "  Starting point of the grid  = " << keys.rad_str << endl;
        cout << "  Grid step parameters        = " << keys.rad_stp << endl;
    };
    cout << endl;
    ifile.close();
    
    line.clear();
    ifile.open ( inpname );
    while( getline(ifile, line) ) if( line == "$CAPAR" ) { 
        getline(ifile, line); keys.cap_on  = atof( line.c_str() );
        getline(ifile, line); keys.cap_mul = atof( line.c_str() );
        break;
    };
    
    if( keys.cap_on  < 0.0 || 
        keys.cap_mul < 0.0  ) keys.capi = false; /// fail-safe
    
    if( keys.capi ) {
        cout << " Complex absorbing potentials parameters:" << endl;
        cout << "  Onset of the potential   = " << keys.cap_on  << endl;
        cout << "  Multiplicative parameter = " << keys.cap_mul << endl;
    };
    cout << endl;
    ifile.close();
    
    /* ---------------------------------- */
    /* ----- ONE-ELECTRON INTEGRALS ----- */
    /* ---------------------------------- */
    
    int num_crt2,num_sph2;
    int shl_siz_crt,shl_siz_crt2,shl_siz_sph,shl_siz_sph2;
    num_crt2 = num_crt * num_crt;
    num_sph2 = num_sph * num_sph;
    
    shl_siz_crt  = crt_siz[bas_lmax];
    shl_siz_sph  = sph_siz[bas_lmax];
    shl_siz_crt2 = shl_siz_crt * shl_siz_crt;
    shl_siz_sph2 = shl_siz_sph * shl_siz_sph;
    
    /* memory for complete matrices */
    arr1d <cdouble> ovrl_crt ( num_crt2 );
    arr1d <cdouble> ovrl_sph ( num_crt2 );
    
    arr1d <cdouble> kin_crt  ( num_crt2 );
    arr1d <cdouble> kin_sph  ( num_crt2 );
    
    arr1d <cdouble> nuc_crt  ( num_crt2 );
    arr1d <cdouble> nuc_sph  ( num_crt2 );
    
    arr1d <cdouble> bare_crt ( num_crt2 );
    arr1d <cdouble> bare_sph ( num_crt2 );
    
    arr1d <cdouble> dipx_crt ( num_crt2 );
    arr1d <cdouble> dipx_sph ( num_crt2 );
    
    arr1d <cdouble> dipy_crt ( num_crt2 );
    arr1d <cdouble> dipy_sph ( num_crt2 );
    
    arr1d <cdouble> dipz_crt ( num_crt2 );
    arr1d <cdouble> dipz_sph ( num_crt2 );
    
    arr1d <cdouble> qdxx_crt ( num_crt2 );
    arr1d <cdouble> qdxx_sph ( num_crt2 );
    
    arr1d <cdouble> qdyy_crt ( num_crt2 );
    arr1d <cdouble> qdyy_sph ( num_crt2 );
    
    arr1d <cdouble> qdzz_crt ( num_crt2 );
    arr1d <cdouble> qdzz_sph ( num_crt2 );
    
    arr1d <cdouble> qdxy_crt ( num_crt2 );
    arr1d <cdouble> qdxy_sph ( num_crt2 );
    
    arr1d <cdouble> qdxz_crt ( num_crt2 );
    arr1d <cdouble> qdxz_sph ( num_crt2 );
    
    arr1d <cdouble> qdyz_crt ( num_crt2 );
    arr1d <cdouble> qdyz_sph ( num_crt2 );
    
    arr1d <cdouble> grdx_crt ( num_crt2 );
    arr1d <cdouble> grdx_sph ( num_crt2 );
    
    arr1d <cdouble> grdy_crt ( num_crt2 );
    arr1d <cdouble> grdy_sph ( num_crt2 );
    
    arr1d <cdouble> grdz_crt ( num_crt2 );
    arr1d <cdouble> grdz_sph ( num_crt2 );
    
    arr1d <cdouble> tx_crt   ( num_crt2 );
    arr1d <cdouble> tx_sph   ( num_crt2 );
    
    arr1d <cdouble> ty_crt   ( num_crt2 );
    arr1d <cdouble> ty_sph   ( num_crt2 );
    
    arr1d <cdouble> tz_crt   ( num_crt2 );
    arr1d <cdouble> tz_sph   ( num_crt2 );
    
    arr1d <cdouble> angz_crt ( num_crt2 );
    arr1d <cdouble> angz_sph ( num_crt2 );
    
    arr1d <cdouble> capi_crt ( num_crt2 );
    arr1d <cdouble> capi_sph ( num_crt2 );
    
    /* memory for shell data */
    arr1d <cdouble> shl_crt_s ( shl_siz_crt2 );
    arr1d <cdouble> shl_sph_s ( shl_siz_sph2 );
    
    arr1d <cdouble> shl_crt_t ( shl_siz_crt2 );
    arr1d <cdouble> shl_sph_t ( shl_siz_sph2 );
    
    arr1d <cdouble> shl_crt_v ( shl_siz_crt2 );
    arr1d <cdouble> shl_sph_v ( shl_siz_sph2 );
    
    arr1d <cdouble> shl_crt_dx( shl_siz_crt2 );
    arr1d <cdouble> shl_sph_dx( shl_siz_sph2 );
    
    arr1d <cdouble> shl_crt_dy( shl_siz_crt2 );
    arr1d <cdouble> shl_sph_dy( shl_siz_sph2 );
    
    arr1d <cdouble> shl_crt_dz( shl_siz_crt2 );
    arr1d <cdouble> shl_sph_dz( shl_siz_sph2 );
    
    arr1d <cdouble> shl_crt_xx( shl_siz_crt2 );
    arr1d <cdouble> shl_sph_xx( shl_siz_sph2 );
    
    arr1d <cdouble> shl_crt_yy( shl_siz_crt2 );
    arr1d <cdouble> shl_sph_yy( shl_siz_sph2 );
    
    arr1d <cdouble> shl_crt_zz( shl_siz_crt2 );
    arr1d <cdouble> shl_sph_zz( shl_siz_sph2 );
    
    arr1d <cdouble> shl_crt_xy( shl_siz_crt2 );
    arr1d <cdouble> shl_sph_xy( shl_siz_sph2 );
    
    arr1d <cdouble> shl_crt_xz( shl_siz_crt2 );
    arr1d <cdouble> shl_sph_xz( shl_siz_sph2 );
    
    arr1d <cdouble> shl_crt_yz( shl_siz_crt2 );
    arr1d <cdouble> shl_sph_yz( shl_siz_sph2 );
    
    arr1d <cdouble> shl_crt_gx( shl_siz_crt2 );
    arr1d <cdouble> shl_sph_gx( shl_siz_sph2 );
    
    arr1d <cdouble> shl_crt_gy( shl_siz_crt2 );
    arr1d <cdouble> shl_sph_gy( shl_siz_sph2 );
    
    arr1d <cdouble> shl_crt_gz( shl_siz_crt2 );
    arr1d <cdouble> shl_sph_gz( shl_siz_sph2 );
    
    arr1d <cdouble> shl_crt_tx( shl_siz_crt2 );
    arr1d <cdouble> shl_sph_tx( shl_siz_sph2 );
    
    arr1d <cdouble> shl_crt_ty( shl_siz_crt2 );
    arr1d <cdouble> shl_sph_ty( shl_siz_sph2 );
    
    arr1d <cdouble> shl_crt_tz( shl_siz_crt2 );
    arr1d <cdouble> shl_sph_tz( shl_siz_sph2 );
    
    arr1d <cdouble> shl_crt_lz( shl_siz_crt2 );
    arr1d <cdouble> shl_sph_lz( shl_siz_sph2 );
    
    arr1d <cdouble> shl_crt_ca( shl_siz_crt2 );
    arr1d <cdouble> shl_sph_ca( shl_siz_sph2 );
    
    arr1d <cdouble> shl_crt_dum ( shl_siz_crt2 );
    arr1d <cdouble> shl_sph_dum ( shl_siz_sph2 );
    
    /* basic stuff for projections, CAPs, etc. */
    int ltot_max  = 30;
    arr2d <double> alpha_ij( ltot_max + 1, ltot_max + 1 );
    BuildAlphaIJ( ltot_max, alpha_ij.v );
    
    /* stuff for various property integrals */
    constexpr int dim_hash = 7;
    cdouble hash_fill[dim_hash][dim_hash] = { 0.0 };
    int hash[dim_hash][3] = { 0 };
    
    hash[0][0] = -1;
    hash[0][1] = +0;
    hash[0][2] = +0;
    
    hash[1][0] = +0;
    hash[1][1] = -1;
    hash[1][2] = +0;
    
    hash[2][0] = +0;
    hash[2][1] = +0;
    hash[2][2] = -1;
    
    hash[3][0] = +1;
    hash[3][1] = +0;
    hash[3][2] = +0;
    
    hash[4][0] = +0;
    hash[4][1] = +1;
    hash[4][2] = +0;
    
    hash[5][0] = +0;
    hash[5][1] = +0;
    hash[5][2] = +1;
    
    hash[6][0] = +0;
    hash[6][1] = +0;
    hash[6][2] = +0;
    
    int li,lj,lij,mi,mj,mjmx;
    int mi_gam,mj_gam;
    int indi,indj,mig,mjg;
    
    int posimx_crt,posjmx_crt;
    int posimx_sph,posjmx_sph;
    int shgi,shgj,shi,shj;
    int ia,ja,ka,ib,jb,kb;
    int ia1,ja1,ka1,ib1,jb1,kb1;
    
    double dum,ai,aj,aij;
    double Ax,Ay,Az,Bx,By,Bz;
    double ABx,ABy,ABz;
    double kPx,kPy,kPz,kP2;
    double Px,Py,Pz,kPP;
    double APx,APy,APz;
    double BPx,BPy,BPz;
    
    cdouble di,dj,dij,II;
    cdouble fac,Ex0,Ey0,Ez0;
    cdouble Ex0_1,Ey0_1,Ez0_1,Exyz;
    cdouble tmp1s,tmp2s,tmp3s,tmp4s;
    
    dum = 0.0;
    
    posimx_crt = 0;
    posimx_sph = 0;
    /* first loop over the shells */
    for(usint i=0; i<basis.size(); i++) {
        li = basis[i].lA;
        Ax = basis[i].Ax;
        Ay = basis[i].Ay;
        Az = basis[i].Az;
        
        shi    = sph_siz[li];
        shgi   = crt_siz[li];
        
        posjmx_crt = 0;
        posjmx_sph = 0;
        /* second loop over the shells */
        for(usint j=0; j<=i; j++) {
            lj = basis[j].lA;
            Bx = basis[j].Ax;
            By = basis[j].Ay;
            Bz = basis[j].Az;
            
            kPx = basis[j].kx - basis[i].kx;
            kPy = basis[j].ky - basis[i].ky;
            kPz = basis[j].kz - basis[i].kz;
            kP2 = kPx * kPx + kPy * kPy + kPz * kPz;
            
            ABx = Ax - Bx;
            ABy = Ay - By;
            ABz = Az - Bz;
            
            /* momentum transfer coefficients */
            Ex0 = exp( I * ( basis[i].kx * Ax - basis[j].kx * Bx ) );
            Ey0 = exp( I * ( basis[i].ky * Ay - basis[j].ky * By ) );
            Ez0 = exp( I * ( basis[i].kz * Az - basis[j].kz * Bz ) );
            
            lij    = li + lj;
            shj    = sph_siz[lj];
            shgj   = crt_siz[lj];
            
            ECoefs <double> Eijx( li + 1, lj + 1, dum, dum, ABx );
            ECoefs <double> Eijy( li + 1, lj + 1, dum, dum, ABy );
            ECoefs <double> Eijz( li + 1, lj + 1, dum, dum, ABz );
            
            RInts1E <cdouble> ROvrl( lij + 2 );
            RInts1E <cdouble> RNucA( lij );
            
            /* CAP factors */
            arr2d < double> cap_ang( lij + 1, ( lij + 1 ) * ( lij + 2 ) / 2 );
            arr2d < double> cap_rad( lij + 1, lij + 5 );
            arr1d <cdouble> cap_tmp( ( lij + 1 ) * ( lij + 2 ) / 2 );
            
            if( keys.capi )
                CAPAngFac( kPx, kPy, kPz, lij, lij, alpha_ij.v, cap_ang.v );
            
            /* zero out the shell data */
            shl_crt_s.zero(); shl_sph_s.zero();
            shl_crt_t.zero(); shl_sph_t.zero();
            shl_crt_v.zero(); shl_sph_v.zero();
            
            shl_crt_dx.zero(); shl_sph_dx.zero();
            shl_crt_dy.zero(); shl_sph_dy.zero();
            shl_crt_dz.zero(); shl_sph_dz.zero();
            
            shl_crt_xx.zero(); shl_sph_xx.zero();
            shl_crt_yy.zero(); shl_sph_yy.zero();
            shl_crt_zz.zero(); shl_sph_zz.zero();
            shl_crt_xy.zero(); shl_sph_xy.zero();
            shl_crt_xz.zero(); shl_sph_xz.zero();
            shl_crt_yz.zero(); shl_sph_yz.zero();
            
            shl_crt_gx.zero(); shl_sph_gx.zero();
            shl_crt_gy.zero(); shl_sph_gy.zero();
            shl_crt_gz.zero(); shl_sph_gz.zero();
            
            shl_crt_tx.zero(); shl_sph_tx.zero();
            shl_crt_ty.zero(); shl_sph_ty.zero();
            shl_crt_tz.zero(); shl_sph_tz.zero();
            
            shl_crt_lz.zero(); shl_sph_lz.zero();
            shl_crt_ca.zero(); shl_sph_ca.zero();
            
            /* loop over contractions - shell I */
            for(usint k1=0; k1<basis[i].alphaA.size(); k1++ ) {
                ai = basis[i].alphaA[k1];
                di = basis[i].dA_re[k1] - I * basis[i].dA_im[k1];
                
                /* loop over contractions - shell J */
                for(usint k2=0; k2<basis[j].alphaA.size(); k2++ ) {
                    aj = basis[j].alphaA[k2];
                    dj = basis[j].dA_re[k2] + I * basis[j].dA_im[k2];
                    
                    aij = ai + aj;
                    dij = di * dj;
                    
                    Px = ai * Ax + aj * Bx; Px = Px / aij;
                    Py = ai * Ay + aj * By; Py = Py / aij;
                    Pz = ai * Az + aj * Bz; Pz = Pz / aij;
                    
                    Eijx.zero();
                    Eijy.zero();
                    Eijz.zero();
                    
                    Eijx.load( ai, aj );
                    Eijy.load( ai, aj );
                    Eijz.load( ai, aj );
                    
                    /* calculate the Eij coefficients */
                    CalcEijt( Eijx );
                    CalcEijt( Eijy );
                    CalcEijt( Eijz );
                    
                    Ex0_1 = Ex0 * exp( - ai * aj * ABx * ABx / aij );
                    Ey0_1 = Ey0 * exp( - ai * aj * ABy * ABy / aij );
                    Ez0_1 = Ez0 * exp( - ai * aj * ABz * ABz / aij );
                    Exyz  = Ex0_1 * Ey0_1 * Ez0_1;
                    
                    kPP = kPx * Px + kPy * Py + kPz * Pz;
                    fac = pow( M_PI / aij, 1.5 ) 
                        * exp( -0.25 * kP2 / aij )
                        * exp( I * kPP );
                    
                    /* calculate R_tuv integrals - overlap/kinetic */
                    ROvrl.zero();
                    ROvrl.load( kPx, kPy, kPz, Px, Py, Pz, aij );
                    CalcROvrl ( ROvrl, fac );
                    
                    /* radial CAP integrals */
                    cap_rad.zero(); cap_tmp.zero();
                    if( keys.capi ) {
                        CAPRadFac( lij, lij + 4, sqrt( kP2 ) * keys.cap_on,
                                   aij * keys.cap_on * keys.cap_on, cap_rad.v );
                        
                        ///for(int l=0; l<=lij; l++) for(int n=l+2; n<=lij+4; n++)
                        ///    cout << l << " " << n << " " << cap_rad.v[l][n] << endl;
                        
                        usint mij;
                        double pre_fac;
                        pre_fac = 4.0 * M_PI * pow( keys.cap_on, lij + 5 );
                        
                        for(mij=0; mij<cap_tmp.num; mij++) {
                            II = 1.0;
                            for(usint L=0; L<=lij; L++) {
                                cap_tmp[mij] += II * cap_ang.v[L][mij] *
                                ( cap_rad.v[L][lij+4] - 2.0 * cap_rad.v[L][lij+3]
                                + cap_rad.v[L][lij+2] );
                                II = II * I;
                            };
                            cap_tmp[mij] = cap_tmp[mij] * pre_fac;
                        };
                    };
                    
                    ia = 0;
                    ja = 0;
                    /* loop over Cartesian components - shell I */
                    for(mi=0; mi<shgi; mi++) {
                        ka = li - ia - ja;
                        ib = 0;
                        jb = 0;
                        /* loop over Cartesian components - shell J */
                        for(mj=0; mj<shgj; mj++) {
                            kb = lj - ib - jb;
                            
                            /* filling the hash table */
                            for(int hi=0; hi<dim_hash; hi++) {
                                ia1 = ia + hash[hi][0];
                                ja1 = ja + hash[hi][1];
                                ka1 = ka + hash[hi][2];
                                if( ia1 < 0 || ja1 < 0 || ka1 < 0 ) continue ;
                                
                                for(int hj=0; hj<dim_hash; hj++) {
                                    ib1 = ib + hash[hj][0];
                                    jb1 = jb + hash[hj][1];
                                    kb1 = kb + hash[hj][2];
                                    if( ib1 < 0 || jb1 < 0 || kb1 < 0 ) continue ;
                                    
                                    tmp3s = 0.0;
                                    for(int t=0; t<=ia1+ib1; t++) {
                                        tmp2s = 0.0;
                                        for(int u=0; u<=ja1+jb1; u++) {
                                            tmp1s = 0.0;
                                            for(int v=0; v<=ka1+kb1; v++) {
                                                tmp1s += ROvrl.rtuv[t][u][v] * Eijz.v[ka1][kb1][v];
                                            };
                                            tmp2s += tmp1s * Eijy.v[ja1][jb1][u];
                                        };
                                        tmp3s += tmp2s * Eijx.v[ia1][ib1][t];
                                    };
                                    
                                    hash_fill[hi][hj] = tmp3s;
                                };
                            };
                            
                            
                            /* overlap integrals */
                            shl_crt_s[ mi * shgj + mj ] += hash_fill[6][6] * dij * Exyz;
                            
                            /* kinetic energy integrals */
                            tmp1s = 0.0; tmp2s = 0.0; tmp3s = 0.0;
                            
                            tmp1s = ( hash_fill[0][0] * ib - 2.0 * aj * hash_fill[0][3] + I * basis[j].kx * hash_fill[0][6] ) * ia
                                  - ( hash_fill[3][0] * ib - 2.0 * aj * hash_fill[3][3] + I * basis[j].kx * hash_fill[3][6] ) * 2.0 * ai
                                  - ( hash_fill[6][0] * ib - 2.0 * aj * hash_fill[6][3] + I * basis[j].kx * hash_fill[6][6] ) * I * basis[i].kx;
                            
                            tmp2s = ( hash_fill[1][1] * jb - 2.0 * aj * hash_fill[1][4] + I * basis[j].ky * hash_fill[1][6] ) * ja
                                  - ( hash_fill[4][1] * jb - 2.0 * aj * hash_fill[4][4] + I * basis[j].ky * hash_fill[4][6] ) * 2.0 * ai
                                  - ( hash_fill[6][1] * jb - 2.0 * aj * hash_fill[6][4] + I * basis[j].ky * hash_fill[6][6] ) * I * basis[i].ky;
                            
                            tmp3s = ( hash_fill[2][2] * kb - 2.0 * aj * hash_fill[2][5] + I * basis[j].kz * hash_fill[2][6] ) * ka
                                  - ( hash_fill[5][2] * kb - 2.0 * aj * hash_fill[5][5] + I * basis[j].kz * hash_fill[5][6] ) * 2.0 * ai
                                  - ( hash_fill[6][2] * kb - 2.0 * aj * hash_fill[6][5] + I * basis[j].kz * hash_fill[6][6] ) * I * basis[i].kz;
                            
                            shl_crt_t [ mi * shgj + mj ] += ( tmp1s + tmp2s + tmp3s ) * dij * Exyz / 2.0;
                            
                            shl_crt_tx[ mi * shgj + mj ] += tmp1s * dij * Exyz / 2.0;
                            shl_crt_ty[ mi * shgj + mj ] += tmp2s * dij * Exyz / 2.0;
                            shl_crt_tz[ mi * shgj + mj ] += tmp3s * dij * Exyz / 2.0;
                            
                            /* dipole moment integrals */
                            if( keys.dip ) {
                                BPx = basis[j].Ax - std::get<0>( keys.points[0] );
                                BPy = basis[j].Ay - std::get<1>( keys.points[0] );
                                BPz = basis[j].Az - std::get<2>( keys.points[0] );
                            
                                tmp1s = hash_fill[6][3] + BPx * hash_fill[6][6];
                                tmp2s = hash_fill[6][4] + BPy * hash_fill[6][6];
                                tmp3s = hash_fill[6][5] + BPz * hash_fill[6][6];
                                
                                shl_crt_dx[ mi * shgj + mj ] += tmp1s * dij * Exyz;
                                shl_crt_dy[ mi * shgj + mj ] += tmp2s * dij * Exyz;
                                shl_crt_dz[ mi * shgj + mj ] += tmp3s * dij * Exyz;
                            };
                            
                            /* quadruple moment integrals */
                            if( keys.quad ) {
                                APx = basis[i].Ax - std::get<0>( keys.points[0] );
                                APy = basis[i].Ay - std::get<1>( keys.points[0] );
                                APz = basis[i].Az - std::get<2>( keys.points[0] );
                                
                                BPx = basis[j].Ax - std::get<0>( keys.points[0] );
                                BPy = basis[j].Ay - std::get<1>( keys.points[0] );
                                BPz = basis[j].Az - std::get<2>( keys.points[0] );
                                
                                /* x-x */
                                tmp1s = hash_fill[3][3] + hash_fill[3][6] * BPx
                                      + hash_fill[6][3] * APx + hash_fill[6][6] * APx * BPx;
                                
                                shl_crt_xx[ mi * shgj + mj ] += tmp1s * dij * Exyz;
                                
                                /* x-y */
                                tmp1s = hash_fill[3][4] + hash_fill[3][6] * BPy
                                      + hash_fill[6][4] * APx + hash_fill[6][6] * APx * BPy;
                                
                                shl_crt_xy[ mi * shgj + mj ] += tmp1s * dij * Exyz;
                                
                                /* x-z */
                                tmp1s = hash_fill[3][5] + hash_fill[3][6] * BPz
                                      + hash_fill[6][5] * APx + hash_fill[6][6] * APx * BPz;
                                
                                shl_crt_xz[ mi * shgj + mj ] += tmp1s * dij * Exyz;
                                
                                /* y-y */
                                tmp1s = hash_fill[4][4] + hash_fill[4][6] * BPy
                                      + hash_fill[6][4] * APy + hash_fill[6][6] * APy * BPy;
                                
                                shl_crt_yy[ mi * shgj + mj ] += tmp1s * dij * Exyz;
                                
                                /* y-z */
                                tmp1s = hash_fill[4][5] + hash_fill[4][6] * BPz
                                      + hash_fill[6][5] * APy + hash_fill[6][6] * APy * BPz;
                                
                                shl_crt_yz[ mi * shgj + mj ] += tmp1s * dij * Exyz;
                                
                                /* z-z */
                                tmp1s = hash_fill[5][5] + hash_fill[5][6] * BPz
                                      + hash_fill[6][5] * APz + hash_fill[6][6] * APz * BPz;
                                
                                shl_crt_zz[ mi * shgj + mj ] += tmp1s * dij * Exyz;
                            };
                            
                            /* dipole velocity integrals */
                            if( keys.velo ) {
                                tmp1s = hash_fill[6][0] * ib - 2.0 * aj * hash_fill[6][3] + I * basis[j].kx * hash_fill[6][6];
                                tmp2s = hash_fill[6][1] * jb - 2.0 * aj * hash_fill[6][4] + I * basis[j].ky * hash_fill[6][6];
                                tmp3s = hash_fill[6][2] * kb - 2.0 * aj * hash_fill[6][5] + I * basis[j].kz * hash_fill[6][6];
                                
                                shl_crt_gx[ mi * shgj + mj ] += tmp1s * dij * Exyz;
                                shl_crt_gy[ mi * shgj + mj ] += tmp2s * dij * Exyz;
                                shl_crt_gz[ mi * shgj + mj ] += tmp3s * dij * Exyz;
                            };
                            
                            /* angular momentum operator, L_z */
                            
                            
                            
                            
                            
                            
                            /* complex absorbing potentials integrals */
                            if( keys.capi ) {
                                usint pos;
                                pos = FindPos( ia + ib, ja + jb, ka + kb );
                                shl_crt_ca[ mi * shgj + mj ] += dij * keys.cap_mul * cap_tmp[pos];
                            };
                            
                            jb++;
                            if( jb > lj - ib ) { jb = 0; ib++; };
                        };
                        ja++;
                        if( ja > li - ia ) { ja = 0; ia++; };
                    };
                    
                    /* loop over all nuclei */
                    for(usint nc = 0; nc<nucli.size(); nc++ ) {
                    
                        /* calculate R_tuv integrals - nuclear attr. */
                        RNucA.zero();
                        RNucA.load( kPx, kPy, kPz, Px, Py, Pz, 
                                    nucli[nc].Cx, nucli[nc].Cy,
                                    nucli[nc].Cz, aij );
                        CalcRNucA ( RNucA );
                        
                        ia = 0;
                        ja = 0;
                        /* loop over Cartesian components - shell I */
                        for(mi=0; mi<shgi; mi++) {
                            ka = li - ia - ja;
                            ib = 0;
                            jb = 0;
                            /* loop over Cartesian components - shell J */
                            for(mj=0; mj<shgj; mj++) {
                                kb = lj - ib - jb;
                                
                                tmp3s = 0.0;
                                for(int t=0; t<=ia+ib; t++) {
                                    tmp2s = 0.0;
                                    for(int u=0; u<=ja+jb; u++) {
                                        tmp1s = 0.0;
                                        for(int v=0; v<=ka+kb; v++) {
                                            tmp1s += RNucA.rtuv[t][u][v] * Eijz.v[ka][kb][v];
                                        };
                                        tmp2s += tmp1s * Eijy.v[ja][jb][u];
                                    };
                                    tmp3s += tmp2s * Eijx.v[ia][ib][t];
                                };
                                
                                shl_crt_v[ mi * shgj + mj ] -= tmp3s * dij * Exyz * nucli[nc].chrg;
                                
                                jb++;
                                if( jb > lj - ib ) { jb = 0; ib++; };
                            };
                            ja++;
                            if( ja > li - ia ) { ja = 0; ia++; };
                        };  
                        
                    }; /* end loop over nuclei */
                };
            }; /* end loop over contractions */
            
            /* rescale to fix the angular normalisation */
            ia = 0;
            ja = 0;
            for(mi=0; mi<shgi; mi++) {
                ka = li - ia - ja;
                ib = 0;
                jb = 0;
                for(mj=0; mj<shgj; mj++) {
                    kb = lj - ib - jb;
                    shl_crt_s[ mi * shgj + mj ] = shl_crt_s[ mi * shgj + mj ] / xyz_norm[li][mi] / xyz_norm[lj][mj];
                    shl_crt_v[ mi * shgj + mj ] = shl_crt_v[ mi * shgj + mj ] / xyz_norm[li][mi] / xyz_norm[lj][mj];
                    shl_crt_t[ mi * shgj + mj ] = shl_crt_t[ mi * shgj + mj ] / xyz_norm[li][mi] / xyz_norm[lj][mj];
                    
                    if( keys.dip ) {
                        shl_crt_dx[ mi * shgj + mj ] = shl_crt_dx[ mi * shgj + mj ] / xyz_norm[li][mi] / xyz_norm[lj][mj];
                        shl_crt_dy[ mi * shgj + mj ] = shl_crt_dy[ mi * shgj + mj ] / xyz_norm[li][mi] / xyz_norm[lj][mj];
                        shl_crt_dz[ mi * shgj + mj ] = shl_crt_dz[ mi * shgj + mj ] / xyz_norm[li][mi] / xyz_norm[lj][mj];
                    };
                    
                    if( keys.quad ) {
                        shl_crt_xx[ mi * shgj + mj ] = shl_crt_xx[ mi * shgj + mj ] / xyz_norm[li][mi] / xyz_norm[lj][mj];
                        shl_crt_yy[ mi * shgj + mj ] = shl_crt_yy[ mi * shgj + mj ] / xyz_norm[li][mi] / xyz_norm[lj][mj];
                        shl_crt_zz[ mi * shgj + mj ] = shl_crt_zz[ mi * shgj + mj ] / xyz_norm[li][mi] / xyz_norm[lj][mj];
                        shl_crt_xy[ mi * shgj + mj ] = shl_crt_xy[ mi * shgj + mj ] / xyz_norm[li][mi] / xyz_norm[lj][mj];
                        shl_crt_xz[ mi * shgj + mj ] = shl_crt_xz[ mi * shgj + mj ] / xyz_norm[li][mi] / xyz_norm[lj][mj];
                        shl_crt_yz[ mi * shgj + mj ] = shl_crt_yz[ mi * shgj + mj ] / xyz_norm[li][mi] / xyz_norm[lj][mj];
                    };
                    
                    if( keys.velo ) {
                        shl_crt_gx[ mi * shgj + mj ] = shl_crt_gx[ mi * shgj + mj ] / xyz_norm[li][mi] / xyz_norm[lj][mj];
                        shl_crt_gy[ mi * shgj + mj ] = shl_crt_gy[ mi * shgj + mj ] / xyz_norm[li][mi] / xyz_norm[lj][mj];
                        shl_crt_gz[ mi * shgj + mj ] = shl_crt_gz[ mi * shgj + mj ] / xyz_norm[li][mi] / xyz_norm[lj][mj];
                    };
                    
                    if( keys.kinc ) {
                        shl_crt_tx[ mi * shgj + mj ] = shl_crt_tx[ mi * shgj + mj ] / xyz_norm[li][mi] / xyz_norm[lj][mj];
                        shl_crt_ty[ mi * shgj + mj ] = shl_crt_ty[ mi * shgj + mj ] / xyz_norm[li][mi] / xyz_norm[lj][mj];
                        shl_crt_tz[ mi * shgj + mj ] = shl_crt_tz[ mi * shgj + mj ] / xyz_norm[li][mi] / xyz_norm[lj][mj];
                    };
                    
                    if( keys.capi )
                        shl_crt_ca[ mi * shgj + mj ] = shl_crt_ca[ mi * shgj + mj ] / xyz_norm[li][mi] / xyz_norm[lj][mj];
                    
                    jb++;
                    if( jb > lj - ib ) { jb = 0; ib++; };
                };
                ja++;
                if( ja > li - ia ) { ja = 0; ia++; };
            };
            
            /* transform to the spherical representation */
            TransToSpher( shl_crt_s.v, shl_sph_s.v, li, lj );
            TransToSpher( shl_crt_t.v, shl_sph_t.v, li, lj );
            TransToSpher( shl_crt_v.v, shl_sph_v.v, li, lj );
            
            if( keys.dip ) {
                TransToSpher( shl_crt_dx.v, shl_sph_dx.v, li, lj );
                TransToSpher( shl_crt_dy.v, shl_sph_dy.v, li, lj );
                TransToSpher( shl_crt_dz.v, shl_sph_dz.v, li, lj );
            };
            
            if( keys.quad ) {
                TransToSpher( shl_crt_xx.v, shl_sph_xx.v, li, lj );
                TransToSpher( shl_crt_yy.v, shl_sph_yy.v, li, lj );
                TransToSpher( shl_crt_zz.v, shl_sph_zz.v, li, lj );
                TransToSpher( shl_crt_xy.v, shl_sph_xy.v, li, lj );
                TransToSpher( shl_crt_xz.v, shl_sph_xz.v, li, lj );
                TransToSpher( shl_crt_yz.v, shl_sph_yz.v, li, lj );
            };
            
            if( keys.velo ) {
                TransToSpher( shl_crt_gx.v, shl_sph_gx.v, li, lj );
                TransToSpher( shl_crt_gy.v, shl_sph_gy.v, li, lj );
                TransToSpher( shl_crt_gz.v, shl_sph_gz.v, li, lj );
            };
            
            if( keys.kinc ) {
                TransToSpher( shl_crt_tx.v, shl_sph_tx.v, li, lj );
                TransToSpher( shl_crt_ty.v, shl_sph_ty.v, li, lj );
                TransToSpher( shl_crt_tz.v, shl_sph_tz.v, li, lj );
            };
            
            if( keys.capi )
                TransToSpher( shl_crt_ca.v, shl_sph_ca.v, li, lj );
            
            /* switch to the Gamess holy order */
            if( gamess_order ) {
            
                shl_crt_dum.zero(); GamessHolyOrder( shl_crt_s.v, shl_crt_dum.v, li, lj );
                memcpy( shl_crt_s.v, shl_crt_dum.v, shl_siz_crt2 * sizeof( cdouble ) );
                
                shl_crt_dum.zero(); GamessHolyOrder( shl_crt_t.v, shl_crt_dum.v, li, lj );
                memcpy( shl_crt_t.v, shl_crt_dum.v, shl_siz_crt2 * sizeof( cdouble ) );
                
                shl_crt_dum.zero(); GamessHolyOrder( shl_crt_v.v, shl_crt_dum.v, li, lj );
                memcpy( shl_crt_v.v, shl_crt_dum.v, shl_siz_crt2 * sizeof( cdouble ) );
                
                if( keys.dip ) {
                    shl_crt_dum.zero(); GamessHolyOrder( shl_crt_dx.v, shl_crt_dum.v, li, lj );
                    memcpy( shl_crt_dx.v, shl_crt_dum.v, shl_siz_crt2 * sizeof( cdouble ) );
                    
                    shl_crt_dum.zero(); GamessHolyOrder( shl_crt_dy.v, shl_crt_dum.v, li, lj );
                    memcpy( shl_crt_dy.v, shl_crt_dum.v, shl_siz_crt2 * sizeof( cdouble ) );
                    
                    shl_crt_dum.zero(); GamessHolyOrder( shl_crt_dz.v, shl_crt_dum.v, li, lj );
                    memcpy( shl_crt_dz.v, shl_crt_dum.v, shl_siz_crt2 * sizeof( cdouble ) );
                };
                
                if( keys.quad ) {
                    shl_crt_dum.zero(); GamessHolyOrder( shl_crt_xx.v, shl_crt_dum.v, li, lj );
                    memcpy( shl_crt_xx.v, shl_crt_dum.v, shl_siz_crt2 * sizeof( cdouble ) );
                    
                    shl_crt_dum.zero(); GamessHolyOrder( shl_crt_yy.v, shl_crt_dum.v, li, lj );
                    memcpy( shl_crt_yy.v, shl_crt_dum.v, shl_siz_crt2 * sizeof( cdouble ) );
                    
                    shl_crt_dum.zero(); GamessHolyOrder( shl_crt_zz.v, shl_crt_dum.v, li, lj );
                    memcpy( shl_crt_zz.v, shl_crt_dum.v, shl_siz_crt2 * sizeof( cdouble ) );
                    
                    shl_crt_dum.zero(); GamessHolyOrder( shl_crt_xy.v, shl_crt_dum.v, li, lj );
                    memcpy( shl_crt_xy.v, shl_crt_dum.v, shl_siz_crt2 * sizeof( cdouble ) );
                    
                    shl_crt_dum.zero(); GamessHolyOrder( shl_crt_xz.v, shl_crt_dum.v, li, lj );
                    memcpy( shl_crt_xz.v, shl_crt_dum.v, shl_siz_crt2 * sizeof( cdouble ) );
                    
                    shl_crt_dum.zero(); GamessHolyOrder( shl_crt_yz.v, shl_crt_dum.v, li, lj );
                    memcpy( shl_crt_yz.v, shl_crt_dum.v, shl_siz_crt2 * sizeof( cdouble ) );
                };
                
                if( keys.velo ) {
                    shl_crt_dum.zero(); GamessHolyOrder( shl_crt_gx.v, shl_crt_dum.v, li, lj );
                    memcpy( shl_crt_gx.v, shl_crt_dum.v, shl_siz_crt2 * sizeof( cdouble ) );
                    
                    shl_crt_dum.zero(); GamessHolyOrder( shl_crt_gy.v, shl_crt_dum.v, li, lj );
                    memcpy( shl_crt_gy.v, shl_crt_dum.v, shl_siz_crt2 * sizeof( cdouble ) );
                    
                    shl_crt_dum.zero(); GamessHolyOrder( shl_crt_gz.v, shl_crt_dum.v, li, lj );
                    memcpy( shl_crt_gz.v, shl_crt_dum.v, shl_siz_crt2 * sizeof( cdouble ) );
                };
                
                if( keys.kinc ) {
                    shl_crt_dum.zero(); GamessHolyOrder( shl_crt_tx.v, shl_crt_dum.v, li, lj );
                    memcpy( shl_crt_tx.v, shl_crt_dum.v, shl_siz_crt2 * sizeof( cdouble ) );
                    
                    shl_crt_dum.zero(); GamessHolyOrder( shl_crt_ty.v, shl_crt_dum.v, li, lj );
                    memcpy( shl_crt_ty.v, shl_crt_dum.v, shl_siz_crt2 * sizeof( cdouble ) );
                    
                    shl_crt_dum.zero(); GamessHolyOrder( shl_crt_tz.v, shl_crt_dum.v, li, lj );
                    memcpy( shl_crt_tz.v, shl_crt_dum.v, shl_siz_crt2 * sizeof( cdouble ) );
                };
                
                if( keys.capi ) {
                    shl_crt_dum.zero(); GamessHolyOrder( shl_crt_ca.v, shl_crt_dum.v, li, lj );
                    memcpy( shl_crt_ca.v, shl_crt_dum.v, shl_siz_crt2 * sizeof( cdouble ) );
                };
            };
            
            /* place in the table - Cartesian rep. */
            for(mi=0; mi<shgi; mi++) {
                indi = posimx_crt + mi;
                mjmx = ( i == j ) ? mi + 1 : shgj;
                for(mj=0; mj<mjmx; mj++) {
                    indj = posjmx_crt + mj;
                    
                    ovrl_crt[ indi * num_crt + indj ] = shl_crt_s[ mi * shgj + mj ];
                    ovrl_crt[ indj * num_crt + indi ] = conj( shl_crt_s[ mi * shgj + mj ] );
                    
                    nuc_crt [ indi * num_crt + indj ] = shl_crt_v[ mi * shgj + mj ];
                    nuc_crt [ indj * num_crt + indi ] = conj( shl_crt_v[ mi * shgj + mj ] );
                    
                    kin_crt [ indi * num_crt + indj ] = shl_crt_t[ mi * shgj + mj ];
                    kin_crt [ indj * num_crt + indi ] = conj( shl_crt_t[ mi * shgj + mj ] );
                    
                    if( keys.dip ) {
                        dipx_crt[ indi * num_crt + indj ] = shl_crt_dx[ mi * shgj + mj ];
                        dipx_crt[ indj * num_crt + indi ] = conj( shl_crt_dx[ mi * shgj + mj ] );
                        
                        dipy_crt[ indi * num_crt + indj ] = shl_crt_dy[ mi * shgj + mj ];
                        dipy_crt[ indj * num_crt + indi ] = conj( shl_crt_dy[ mi * shgj + mj ] );
                        
                        dipz_crt[ indi * num_crt + indj ] = shl_crt_dz[ mi * shgj + mj ];
                        dipz_crt[ indj * num_crt + indi ] = conj( shl_crt_dz[ mi * shgj + mj ] );
                    };
                    
                    if( keys.quad ) {
                        qdxx_crt[ indi * num_crt + indj ] = shl_crt_xx[ mi * shgj + mj ];
                        qdxx_crt[ indj * num_crt + indi ] = conj( shl_crt_xx[ mi * shgj + mj ] );
                        
                        qdyy_crt[ indi * num_crt + indj ] = shl_crt_yy[ mi * shgj + mj ];
                        qdyy_crt[ indj * num_crt + indi ] = conj( shl_crt_yy[ mi * shgj + mj ] );
                        
                        qdzz_crt[ indi * num_crt + indj ] = shl_crt_zz[ mi * shgj + mj ];
                        qdzz_crt[ indj * num_crt + indi ] = conj( shl_crt_zz[ mi * shgj + mj ] );
                        
                        qdxy_crt[ indi * num_crt + indj ] = shl_crt_xy[ mi * shgj + mj ];
                        qdxy_crt[ indj * num_crt + indi ] = conj( shl_crt_xy[ mi * shgj + mj ] );
                        
                        qdxz_crt[ indi * num_crt + indj ] = shl_crt_xz[ mi * shgj + mj ];
                        qdxz_crt[ indj * num_crt + indi ] = conj( shl_crt_xz[ mi * shgj + mj ] );
                        
                        qdyz_crt[ indi * num_crt + indj ] = shl_crt_yz[ mi * shgj + mj ];
                        qdyz_crt[ indj * num_crt + indi ] = conj( shl_crt_yz[ mi * shgj + mj ] );
                    };
                    
                    if( keys.velo ) {
                        grdx_crt[ indi * num_crt + indj ] = shl_crt_gx[ mi * shgj + mj ];
                        grdx_crt[ indj * num_crt + indi ] = -conj( shl_crt_gx[ mi * shgj + mj ] );
                        
                        grdy_crt[ indi * num_crt + indj ] = shl_crt_gy[ mi * shgj + mj ];
                        grdy_crt[ indj * num_crt + indi ] = -conj( shl_crt_gy[ mi * shgj + mj ] );
                        
                        grdz_crt[ indi * num_crt + indj ] = shl_crt_gz[ mi * shgj + mj ];
                        grdz_crt[ indj * num_crt + indi ] = -conj( shl_crt_gz[ mi * shgj + mj ] );
                    };
                    
                    if( keys.kinc ) {
                        tx_crt  [ indi * num_crt + indj ] = shl_crt_tx[ mi * shgj + mj ];
                        tx_crt  [ indj * num_crt + indi ] = conj( shl_crt_tx[ mi * shgj + mj ] );
                        
                        ty_crt  [ indi * num_crt + indj ] = shl_crt_ty[ mi * shgj + mj ];
                        ty_crt  [ indj * num_crt + indi ] = conj( shl_crt_ty[ mi * shgj + mj ] );
                        
                        tz_crt  [ indi * num_crt + indj ] = shl_crt_tz[ mi * shgj + mj ];
                        tz_crt  [ indj * num_crt + indi ] = conj( shl_crt_tz[ mi * shgj + mj ] );
                    };
                    
                    if( keys.capi ) {
                        capi_crt[ indi * num_crt + indj ] = shl_crt_ca[ mi * shgj + mj ];
                        capi_crt[ indj * num_crt + indi ] = conj( shl_crt_ca[ mi * shgj + mj ] );
                    };
                };
            };
            
            /* place in the table - spherical rep. */
            for(mi=0; mi<shi; mi++) {
                indi = posimx_sph + mi;
                mjmx = ( i == j ) ? mi + 1 : shj;
                for(mj=0; mj<mjmx; mj++) {
                    indj = posjmx_sph + mj;
                    
                    ovrl_sph[ indi * num_sph + indj ] = shl_sph_s[ mi * shj + mj ];
                    ovrl_sph[ indj * num_sph + indi ] = conj( shl_sph_s[ mi * shj + mj ] );
                    
                    nuc_sph [ indi * num_sph + indj ] = shl_sph_v[ mi * shj + mj ];
                    nuc_sph [ indj * num_sph + indi ] = conj( shl_sph_v[ mi * shj + mj ] );
                    
                    kin_sph [ indi * num_sph + indj ] = shl_sph_t[ mi * shj + mj ];
                    kin_sph [ indj * num_sph + indi ] = conj( shl_sph_t[ mi * shj + mj ] );
                    
                    if( keys.dip ) {
                        dipx_sph[ indi * num_sph + indj ] = shl_sph_dx[ mi * shj + mj ];
                        dipx_sph[ indj * num_sph + indi ] = conj( shl_sph_dx[ mi * shj + mj ] );
                        
                        dipy_sph[ indi * num_sph + indj ] = shl_sph_dy[ mi * shj + mj ];
                        dipy_sph[ indj * num_sph + indi ] = conj( shl_sph_dy[ mi * shj + mj ] );
                        
                        dipz_sph[ indi * num_sph + indj ] = shl_sph_dz[ mi * shj + mj ];
                        dipz_sph[ indj * num_sph + indi ] = conj( shl_sph_dz[ mi * shj + mj ] );
                    };
                    
                    if( keys.quad ) {
                        qdxx_sph[ indi * num_sph + indj ] = shl_sph_xx[ mi * shj + mj ];
                        qdxx_sph[ indj * num_sph + indi ] = conj( shl_sph_xx[ mi * shj + mj ] );
                        
                        qdyy_sph[ indi * num_sph + indj ] = shl_sph_yy[ mi * shj + mj ];
                        qdyy_sph[ indj * num_sph + indi ] = conj( shl_sph_yy[ mi * shj + mj ] );
                        
                        qdzz_sph[ indi * num_sph + indj ] = shl_sph_zz[ mi * shj + mj ];
                        qdzz_sph[ indj * num_sph + indi ] = conj( shl_sph_zz[ mi * shj + mj ] );
                        
                        qdxy_sph[ indi * num_sph + indj ] = shl_sph_xy[ mi * shj + mj ];
                        qdxy_sph[ indj * num_sph + indi ] = conj( shl_sph_xy[ mi * shj + mj ] );
                        
                        qdxz_sph[ indi * num_sph + indj ] = shl_sph_xz[ mi * shj + mj ];
                        qdxz_sph[ indj * num_sph + indi ] = conj( shl_sph_xz[ mi * shj + mj ] );
                        
                        qdyz_sph[ indi * num_sph + indj ] = shl_sph_yz[ mi * shj + mj ];
                        qdyz_sph[ indj * num_sph + indi ] = conj( shl_sph_yz[ mi * shj + mj ] );
                    };
                    
                    if( keys.velo ) {
                        grdx_sph[ indi * num_sph + indj ] = shl_sph_gx[ mi * shj + mj ];
                        grdx_sph[ indj * num_sph + indi ] = -conj( shl_sph_gx[ mi * shj + mj ] );
                        
                        grdy_sph[ indi * num_sph + indj ] = shl_sph_gy[ mi * shj + mj ];
                        grdy_sph[ indj * num_sph + indi ] = -conj( shl_sph_gy[ mi * shj + mj ] );
                        
                        grdz_sph[ indi * num_sph + indj ] = shl_sph_gz[ mi * shj + mj ];
                        grdz_sph[ indj * num_sph + indi ] = -conj( shl_sph_gz[ mi * shj + mj ] );
                    };
                    
                    if( keys.kinc ) {
                        tx_sph  [ indi * num_sph + indj ] = shl_sph_tx[ mi * shj + mj ];
                        tx_sph  [ indj * num_sph + indi ] = conj( shl_sph_tx[ mi * shj + mj ] );
                        
                        ty_sph  [ indi * num_sph + indj ] = shl_sph_ty[ mi * shj + mj ];
                        ty_sph  [ indj * num_sph + indi ] = conj( shl_sph_ty[ mi * shj + mj ] );
                        
                        tz_sph  [ indi * num_sph + indj ] = shl_sph_tz[ mi * shj + mj ];
                        tz_sph  [ indj * num_sph + indi ] = conj( shl_sph_tz[ mi * shj + mj ] );
                    };
                    
                    if( keys.capi ) {
                        capi_sph[ indi * num_sph + indj ] = shl_sph_ca[ mi * shj + mj ];
                        capi_sph[ indj * num_sph + indi ] = conj( shl_sph_ca[ mi * shj + mj ] );
                    };
                };
            };
            
            posjmx_crt += shgj;
            posjmx_sph += shj;
        };
        posimx_crt += shgi;
        posimx_sph += shi;
    };
    
    for(int i=0; i<num_crt2; i++) bare_crt[i] = kin_crt[i] + nuc_crt[i];
    for(int i=0; i<num_sph2; i++) bare_sph[i] = kin_sph[i] + nuc_sph[i];
    
    ///PrintCMatrix( nuc_crt.v , num_crt );
    ///PrintCMatrix( ovrl_crt.v, num_crt );
    ///PrintCMatrix( kin_crt.v , num_crt );
    ///PrintCMatrix( bare_crt.v, num_crt );
    ///PrintCMatrix( dipz_crt.v, num_crt );
    ///PrintCMatrix( ovrl_crt.v, num_crt );
    ///PrintCMatrix( capi_crt.v, num_crt );
    ///PrintCMatrix( ovrl_crt.v, num_crt );
    
    /* write the one-electron integrals to the disk - cartesian */
    std::ofstream ofs( keys.file1E, std::ios::out|std::ios::binary );
    
    WriteDown( ovrl_crt.v, num_crt2, ofs );
    WriteDown( kin_crt.v , num_crt2, ofs );
    WriteDown( nuc_crt.v , num_crt2, ofs );
    WriteDown( bare_crt.v, num_crt2, ofs );
    
    WriteDown( dipx_crt.v, num_crt2, ofs );
    WriteDown( dipy_crt.v, num_crt2, ofs );
    WriteDown( dipz_crt.v, num_crt2, ofs );
    
    WriteDown( qdxx_crt.v, num_crt2, ofs );
    WriteDown( qdyy_crt.v, num_crt2, ofs );
    WriteDown( qdzz_crt.v, num_crt2, ofs );
    WriteDown( qdxy_crt.v, num_crt2, ofs );
    WriteDown( qdxz_crt.v, num_crt2, ofs );
    WriteDown( qdyz_crt.v, num_crt2, ofs );
    
    WriteDown( grdx_crt.v, num_crt2, ofs );
    WriteDown( grdy_crt.v, num_crt2, ofs );
    WriteDown( grdz_crt.v, num_crt2, ofs );
    
    WriteDown( tx_crt.v  , num_crt2, ofs );
    WriteDown( ty_crt.v  , num_crt2, ofs );
    WriteDown( tz_crt.v  , num_crt2, ofs );
    
    WriteDown( capi_crt.v, num_crt2, ofs );
    
    ofs.close();
    
    /* write the one-electron integrals to the disk - spherical */
    /*
>>>>>>> master
    std::ofstream ofs( keys.file1E, std::ios::out|std::ios::binary );
    
    WriteDown( ovrl_sph.v, num_sph2, ofs );
    WriteDown( kin_sph.v , num_sph2, ofs );
    WriteDown( nuc_sph.v , num_sph2, ofs );
    WriteDown( bare_sph.v, num_sph2, ofs );
    
    WriteDown( dipx_sph.v, num_sph2, ofs );
    WriteDown( dipy_sph.v, num_sph2, ofs );
    WriteDown( dipz_sph.v, num_sph2, ofs );
     
    WriteDown( qdxx_sph.v, num_sph2, ofs );
    WriteDown( qdyy_sph.v, num_sph2, ofs );
    WriteDown( qdzz_sph.v, num_sph2, ofs );
    WriteDown( qdxy_sph.v, num_sph2, ofs );
    WriteDown( qdxz_sph.v, num_sph2, ofs );
    WriteDown( qdyz_sph.v, num_sph2, ofs );
    
    WriteDown( grdx_sph.v, num_sph2, ofs );
    WriteDown( grdy_sph.v, num_sph2, ofs );
    WriteDown( grdz_sph.v, num_sph2, ofs );
    
    WriteDown( tx_sph.v  , num_sph2, ofs );
    WriteDown( ty_sph.v  , num_sph2, ofs );
    WriteDown( tz_sph.v  , num_sph2, ofs );
    
    WriteDown( capi_sph.v, num_crt2, ofs );
    
    ofs.close();
    */
<<<<<<< HEAD

	// PrintCMatrix( ovrl_crt.v , num_crt );

	cout << endl;
	cout << " All one-electron integrals done." << endl;

	/* projection of the basis set functions */
	if (keys.proj) {
		double rr, rr2, kr, valr, vali;
		cdouble g_rad;

		arr1d<cdouble> sum_in_crt(crt_siz[bas_lmax]);
		arr1d<cdouble> sum_in_sph(sph_siz[bas_lmax]);

		arr2d<cdouble> proj_crt(num_crt, keys.n_rad);
		arr2d<cdouble> proj_sph(num_sph, keys.n_rad);

		arr1d<double> ang_fac(crt_siz[bas_lmax]);
		arr1d<double> jn(ltot_max + 1);

		/* loop over projections first */  /// suboptimal but avoids huge storage
		for (usint lP = 0; lP <= keys.lproj_max; lP++) {
			proj_crt.zero();

			/* loop over radial points */
			for (usint nr = 0; nr < keys.n_rad; nr++) {
				rr = keys.rad_str + keys.rad_stp * nr;
				rr2 = rr * rr;

				posimx_crt = 0;
				posimx_sph = 0;
				/* loop over the shells */
				for (usint i = 0; i < basis.size(); i++) {
					li = basis[i].lA;
					shi = sph_siz[li];
					shgi = crt_siz[li];

					kPx = basis[i].kx;
					kPy = basis[i].ky;
					kPz = basis[i].kz;

					kP2 = kPx * kPx + kPy * kPy + kPz * kPz;
					kP2 = sqrt(kP2);
					kr = kP2 * rr;

					/* loop over contractions */
					g_rad = 0.0;
					for (usint k1 = 0; k1 < basis[i].alphaA.size(); k1++) {
						ai = basis[i].alphaA[k1];
						di = basis[i].dA_re[k1] + I * basis[i].dA_im[k1];
						g_rad += exp(-ai * rr2) * di;
					};
					g_rad = g_rad * fourpi * pow(rr, li);

					/* inner summation over L */
					jn.zero();
					sum_in_crt.zero();
					SphericalBesselJ(lP + li, kr, jn.v);

					II = 1.0;
					for (usint L = 0; L <= lP + li; L++) {
						CompAngFac(kPx, kPy, kPz, lP, li,
						           L, alpha_ij.v, ang_fac.v);
						for (mi = 0; mi < shgi; mi++)
							sum_in_crt[mi] += II * jn.v[L] * ang_fac.v[mi];
						II = II * I;
					};

					/* add the Gaussian factor, fix the norm */
					ia = 0;
					ja = 0;
					for (mi = 0; mi < shgi; mi++) {
						ka = li - ia - ja;
						sum_in_crt[mi] = sum_in_crt[mi] * g_rad / xyz_norm[li][mi];

						ja++;
						if (ja > li - ia) {
							ja = 0;
							ia++;
						};
					};

					/* spherical representation */
					sum_in_sph.zero();
					TransToSpher1D(sum_in_crt.v, sum_in_sph.v, li);

					/* switch to Gamess order, finalise */
					for (mi = 0; mi < shgi; mi++) {
						indi = posimx_crt;
						if (gamess_order)
							indi += pos_change_gamess(li, mi);
						else
							indi += mi;
						proj_crt.v[indi][nr] = sum_in_crt[mi];
					};

					for (mi = 0; mi < shi; mi++) {
						indi = posimx_sph + mi;
						proj_sph.v[indi][nr] = sum_in_sph[mi];
					};

					posimx_crt += shgi;
					posimx_sph += shi;
				};
			};
			/* write down */
			string proj_name;
			proj_name = keys.proj1E + "_L" + to_string(lP) + ".F";
			ofstream ofs_proj(proj_name, std::ios::out | std::ios::binary);
			if (!ofs_proj.is_open()) {
				cout << "Cannot open proj1E!\n";
				return EXIT_FAILURE;
			}

			for (mi = 0; mi < num_crt; mi++)
				for (usint nr = 0; nr < keys.n_rad; nr++) {
					valr = proj_crt.v[mi][nr].real();
					vali = proj_crt.v[mi][nr].imag();

					ofs_proj.write(reinterpret_cast<char *>(&valr), sizeof(double));
					ofs_proj.write(reinterpret_cast<char *>(&vali), sizeof(double));
				};

			/*
=======
    
    cout << endl;
    cout << " All one-electron integrals done." << endl;
    
    /* projection of the basis set functions */
    if( keys.proj ) {
    
    double rr,rr2,kr,valr,vali;
    cdouble g_rad;
    
    arr1d <cdouble> sum_in_crt( crt_siz[bas_lmax] );
    arr1d <cdouble> sum_in_sph( sph_siz[bas_lmax] );
    
    arr2d <cdouble> proj_crt( num_crt, keys.n_rad );
    arr2d <cdouble> proj_sph( num_sph, keys.n_rad );
    
    arr1d < double> ang_fac ( crt_siz[bas_lmax] );
    arr1d < double> jn      ( ltot_max + 1 );
    
    /* loop over projections first */ /// suboptimal but avoids huge storage
    for(usint lP=0; lP<=keys.lproj_max; lP++) { 
        proj_crt.zero();
        
        /* loop over radial points */
        for(usint nr=0; nr<keys.n_rad; nr++) {
            rr  = keys.rad_str + keys.rad_stp * nr;
            rr2 = rr * rr;
            
            posimx_crt = 0;
            posimx_sph = 0;
            /* loop over the shells */
            for(usint i=0; i<basis.size(); i++) {
                li   = basis[i].lA;
                shi  = sph_siz[li];
                shgi = crt_siz[li];
            
                kPx = basis[i].kx;
                kPy = basis[i].ky;
                kPz = basis[i].kz;
                
                kP2 = kPx * kPx + kPy * kPy + kPz * kPz;
                kP2 = sqrt( kP2 );
                kr  = kP2 * rr;
            
                /* loop over contractions */
                g_rad = 0.0;
                for(usint k1=0; k1<basis[i].alphaA.size(); k1++ ) {
                    ai = basis[i].alphaA[k1];
                    di = basis[i].dA_re[k1] + I * basis[i].dA_im[k1];
                    g_rad += exp( -ai * rr2 ) * di;
                };
                g_rad = g_rad * fourpi * pow( rr, li );
                
                /* inner summation over L */
                jn.zero(); sum_in_crt.zero();
                SphericalBesselJ( lP+li, kr, jn.v );
                
                II = 1.0;
                for(usint L=0; L<=lP+li; L++) {
                    CompAngFac( kPx, kPy, kPz, lP, li, 
                                L, alpha_ij.v, ang_fac.v );
                    for(mi=0; mi<shgi; mi++)
                        sum_in_crt[mi] += II * jn.v[L] * ang_fac.v[mi];
                    II = II * I;
                };
                
                /* add the Gaussian factor, fix the norm */
                ia = 0;
                ja = 0;
                for(mi=0; mi<shgi; mi++) {
                    ka = li - ia - ja;
                    sum_in_crt[mi] = sum_in_crt[mi] * g_rad / xyz_norm[li][mi];
                    
                    ja++;
                    if( ja > li - ia ) { ja = 0; ia++; };
                };
                
                /* spherical representation */
                sum_in_sph.zero();
                TransToSpher1D( sum_in_crt.v, sum_in_sph.v, li );
                
                /* switch to Gamess order, finalise */
                for(mi=0; mi<shgi; mi++) {
                    indi = posimx_crt;
                    if( gamess_order ) indi += pos_change_gamess( li, mi );
                        else indi += mi;
                    proj_crt.v[indi][nr] = sum_in_crt[mi];
                };
                
                for(mi=0; mi<shi; mi++) {
                    indi = posimx_sph + mi;
                    proj_sph.v[indi][nr] = sum_in_sph[mi];
                };
                
                posimx_crt += shgi;
                posimx_sph += shi;
            };
        };
        /* write down */
        string proj_name;
        proj_name = keys.proj1E + "_L" + to_string( lP ) + ".F";
        ofstream ofs_proj( proj_name, std::ios::out|std::ios::binary );
        
        for(mi=0; mi<num_crt; mi++) for(usint nr=0; nr<keys.n_rad; nr++) {
            valr = proj_crt.v[mi][nr].real();
            vali = proj_crt.v[mi][nr].imag();
            
            ofs_proj.write(reinterpret_cast<char*>(&valr),sizeof(double));
            ofs_proj.write(reinterpret_cast<char*>(&vali),sizeof(double));
        };
        
        /*
>>>>>>> master
        for(mi=0; mi<num_sph; mi++) for(usint nr=0; nr<keys.n_rad; nr++) {
            valr = proj_sph.v[mi][nr].real();
            vali = proj_sph.v[mi][nr].imag();
            
            ofs_proj.write(reinterpret_cast<char*>(&valr),sizeof(double));
            ofs_proj.write(reinterpret_cast<char*>(&vali),sizeof(double));
        };
        */

			ofs_proj.close();
		};

		cout << " Projection integrals done." << endl;
		cout << endl;
	};
	/* ---------------------------------- */
	/* ----- TWO-ELECTRON INTEGRALS ----- */
	/* ---------------------------------- */
	if (keys.eri) {
		cout << " Two-electron integrals - engaging the primary loop." << endl;

		/* open the two-electron integral file */
		std::ofstream ofs_2E;
		ofs_2E.open(keys.file2E, std::ios::out | std::ios::binary);
		/*check the file1E path (MS)*/
		if (!ofs_2E.is_open()) {
			cout << "Cannot open file2E!\n";
			return EXIT_FAILURE;
		}

		llint tot_2E = 0;

		posimx_crt = 0;
		vector<int> posimx_crt_v(basis.size());

		for (int i = 0; i < basis.size(); ++i) {
			posimx_crt_v[i] = posimx_crt;
			posimx_crt += crt_siz[basis[i].lA];
		}

		#pragma omp parallel shared(keys, tot_2E, ofs_2E, basis)
		{
			#pragma omp single
			{
				cout << " Using " << omp_get_num_threads() << " threads.\n";
			}

			/* loop over the shells I, J, K, L*/
			#pragma omp for schedule(dynamic, 1) collapse(4)
			for (int i = 0; i < basis.size(); i++)
				for (int j = 0; j < basis.size(); j++)
					for (int k = 0; k < basis.size(); k++)
						for (int l = 0; l < basis.size(); l++) {

							if (j > i) continue;
							if (k > i) continue;
							int lmax = (i == k) ? j + 1 : basis.size();
							if (l >= lmax) continue;

							/* shell I */
							auto li = basis[i].lA;
							auto Ax = basis[i].Ax;
							auto Ay = basis[i].Ay;
							auto Az = basis[i].Az;

							auto shgi = crt_siz[li];

							/* shell J */
							auto lj = basis[j].lA;
							auto Bx = basis[j].Ax;
							auto By = basis[j].Ay;
							auto Bz = basis[j].Az;

							auto kPx = basis[j].kx - basis[i].kx;
							auto kPy = basis[j].ky - basis[i].ky;
							auto kPz = basis[j].kz - basis[i].kz;

							auto ABx = Ax - Bx;
							auto ABy = Ay - By;
							auto ABz = Az - Bz;

							auto EABx0 = exp(I * (basis[i].kx * Ax - basis[j].kx * Bx));
							auto EABy0 = exp(I * (basis[i].ky * Ay - basis[j].ky * By));
							auto EABz0 = exp(I * (basis[i].kz * Az - basis[j].kz * Bz));

							auto lij = li + lj;
							auto shgj = crt_siz[lj];

							ECoefs<double> Eijx(li, lj, dum, dum, ABx);
							ECoefs<double> Eijy(li, lj, dum, dum, ABy);
							ECoefs<double> Eijz(li, lj, dum, dum, ABz);

							/* shell K */
							auto lk = basis[k].lA;
							auto Cx = basis[k].Ax;
							auto Cy = basis[k].Ay;
							auto Cz = basis[k].Az;

							auto shgk = crt_siz[lk];

							/* shell L */
							auto ll = basis[l].lA;
							auto Dx = basis[l].Ax;
							auto Dy = basis[l].Ay;
							auto Dz = basis[l].Az;

							auto kQx = basis[l].kx - basis[k].kx;
							auto kQy = basis[l].ky - basis[k].ky;
							auto kQz = basis[l].kz - basis[k].kz;

							auto CDx = Cx - Dx;
							auto CDy = Cy - Dy;
							auto CDz = Cz - Dz;

							auto ECDx0 = exp(I * (basis[i].kx * Ax - basis[j].kx * Bx));
							auto ECDy0 = exp(I * (basis[i].ky * Ay - basis[j].ky * By));
							auto ECDz0 = exp(I * (basis[i].kz * Az - basis[j].kz * Bz));

							auto lkl = lk + ll;
							auto shgl = crt_siz[ll];

							ECoefs<double> Eklx(lk, ll, dum, dum, CDx);
							ECoefs<double> Ekly(lk, ll, dum, dum, CDy);
							ECoefs<double> Eklz(lk, ll, dum, dum, CDz);

							RInts2E<cdouble> R_ERI(lij, lkl);

							/* intermediate storage: half- and fully-contracted */
							HalfC<cdouble> half_trans(lij, lk, ll);
							FullC<cdouble> full_trans_crt(li, lj, lk, ll);
							full_trans_crt.zero();

							/* loop over contractions - shell I */
							for (usint k1 = 0; k1 < basis[i].alphaA.size(); k1++) {
								auto ai = basis[i].alphaA[k1];
								auto di = basis[i].dA_re[k1] - I * basis[i].dA_im[k1];

								/* loop over contractions - shell J */
								for (usint k2 = 0; k2 < basis[j].alphaA.size(); k2++) {
									auto aj = basis[j].alphaA[k2];
									auto dj = basis[j].dA_re[k2] + I * basis[j].dA_im[k2];

									auto aij = ai + aj;
									auto dij = di * dj;

									auto Px = ai * Ax + aj * Bx;
									Px = Px / aij;
									auto Py = ai * Ay + aj * By;
									Py = Py / aij;
									auto Pz = ai * Az + aj * Bz;
									Pz = Pz / aij;

									Eijx.zero();
									Eijy.zero();
									Eijz.zero();

									Eijx.load(ai, aj);
									Eijy.load(ai, aj);
									Eijz.load(ai, aj);

									/* calculate the Eij coefficients */
									CalcEijt(Eijx);
									CalcEijt(Eijy);
									CalcEijt(Eijz);

									auto EABx0_1 = EABx0 * exp(-ai * aj * ABx * ABx / aij);
									auto EABy0_1 = EABy0 * exp(-ai * aj * ABy * ABy / aij);
									auto EABz0_1 = EABz0 * exp(-ai * aj * ABz * ABz / aij);
									auto EABxyz = EABx0_1 * EABy0_1 * EABz0_1;

									half_trans.zero();

									/* loop over contractions - shell K */
									for (usint k3 = 0; k3 < basis[k].alphaA.size(); k3++) {
										auto ak = basis[k].alphaA[k3];
										auto dk = basis[k].dA_re[k3] - I * basis[k].dA_im[k3];

										/* loop over contractions - shell L */
										for (usint k4 = 0; k4 < basis[l].alphaA.size(); k4++) {
											auto al = basis[l].alphaA[k4];
											auto dl = basis[l].dA_re[k4] + I * basis[l].dA_im[k4];

											auto akl = ak + al;
											auto dkl = dk * dl;

											auto Qx = ak * Cx + al * Dx;
											Qx = Qx / akl;
											auto Qy = ak * Cy + al * Dy;
											Qy = Qy / akl;
											auto Qz = ak * Cz + al * Dz;
											Qz = Qz / akl;

											Eklx.zero();
											Ekly.zero();
											Eklz.zero();

											Eklx.load(ak, al);
											Ekly.load(ak, al);
											Eklz.load(ak, al);

											/* calculate the Eij coefficients */
											CalcEijt(Eklx);
											CalcEijt(Ekly);
											CalcEijt(Eklz);

											auto ECDx0_1 = ECDx0 * exp(-ak * al * CDx * CDx / akl);
											auto ECDy0_1 = ECDy0 * exp(-ak * al * CDy * CDy / akl);
											auto ECDz0_1 = ECDz0 * exp(-ak * al * CDz * CDz / akl);
											auto ECDxyz = ECDx0_1 * ECDy0_1 * ECDz0_1;

											/* calculate the R integrals for ERI */
											R_ERI.zero();
											R_ERI.load(kPx, kPy, kPz, kQx, kQy, kQz,
											           Px, Py, Pz, Qx, Qy, Qz, aij, akl);
											CalcRERI(R_ERI);
											///R_ERI.print();

											int ic = 0;
											int jc = 0;
											/* loop over Cartesian components - shell K */
											for (usint mk = 0; mk < shgk; mk++) {
												auto kc = lk - ic - jc;
												int id = 0;
												int jd = 0;
												/* loop over Cartesian components - shell L */
												for (usint ml = 0; ml < shgl; ml++) {
													auto kd = ll - id - jd;

													for (int t1 = 0; t1 <= lij; t1++) {
														for (int u1 = 0; u1 <= lij - t1; u1++) {
															for (int v1 = 0; v1 <= lij - t1 - u1; v1++) {
																cdouble tmp3s = 0.0;
																for (int t2 = 0; t2 <= ic + id; t2++) {
																	cdouble tmp2s = 0.0;
																	for (int u2 = 0; u2 <= jc + jd; u2++) {
																		cdouble tmp1s = 0.0;
																		for (int v2 = 0; v2 <= kc + kd; v2++) {
																			tmp1s += R_ERI.R2E[t1][u1][v1][t2][u2][v2] * Eklz.v[kc][kd][v2];
																		};
																		tmp2s += tmp1s * Ekly.v[jc][jd][u2];
																	};
																	tmp3s += tmp2s * Eklx.v[ic][id][t2];
																};

																half_trans.v[t1][u1][v1][mk][ml] += tmp3s * dkl * ECDxyz;
															};
														};
													};

													jd++;
													if (jd > ll - id) {
														jd = 0;
														id++;
													};
												};
												jc++;
												if (jc > lk - ic) {
													jc = 0;
													ic++;
												};
											};
										};
									};

									int ia = 0;
									int ja = 0;
									/* loop over Cartesian components - shell I */
									for (usint mi = 0; mi < shgi; mi++) {
										auto ka = li - ia - ja;

										int ib = 0;
										int jb = 0;
										/* loop over Cartesian components - shell J */
										for (usint mj = 0; mj < shgj; mj++) {
											auto kb = lj - ib - jb;

											int ic = 0;
											int jc = 0;
											/* loop over Cartesian components - shell K */
											for (usint mk = 0; mk < shgk; mk++) {
												int id = 0;
												int jd = 0;
												/* loop over Cartesian components - shell L */
												for (usint ml = 0; ml < shgl; ml++) {
													cdouble tmp3s = 0.0;
													for (int t1 = 0; t1 <= ia + ib; t1++) {
														cdouble tmp2s = 0.0;
														for (int u1 = 0; u1 <= ja + jb; u1++) {
															cdouble tmp1s = 0.0;
															for (int v1 = 0; v1 <= ka + kb; v1++) {
																tmp1s += half_trans.v[t1][u1][v1][mk][ml] * Eijz.v[ka][kb][v1];
															};
															tmp2s += tmp1s * Eijy.v[ja][jb][u1];
														};
														tmp3s += tmp2s * Eijx.v[ia][ib][t1];
													};

													full_trans_crt.v[mi][mj][mk][ml] += tmp3s * dij * EABxyz;

													jd++;
													if (jd > ll - id) {
														jd = 0;
														id++;
													};
												};
												jc++;
												if (jc > lk - ic) {
													jc = 0;
													ic++;
												};
											};

											jb++;
											if (jb > lj - ib) {
												jb = 0;
												ib++;
											};
										};
										ja++;
										if (ja > li - ia) {
											ja = 0;
											ia++;
										};
									};
								};
							};

							/* rescale to fix the angular normalisation */
							int ia = 0;
							int ja = 0;
							/* loop over Cartesian components - shell I */
							for (usint mi = 0; mi < shgi; mi++) {
								int ib = 0;
								int jb = 0;
								/* loop over Cartesian components - shell J */
								for (usint mj = 0; mj < shgj; mj++) {
									int ic = 0;
									int jc = 0;
									/* loop over Cartesian components - shell K */
									for (usint mk = 0; mk < shgk; mk++) {
										int id = 0;
										int jd = 0;
										/* loop over Cartesian components - shell L */
										for (usint ml = 0; ml < shgl; ml++) {
											full_trans_crt.v[mi][mj][mk][ml] = full_trans_crt.v[mi][mj][mk][ml] / xyz_norm[li][mi] / xyz_norm[lj][mj] / xyz_norm[lk][mk] / xyz_norm[ll][ml];

											jd++;
											if (jd > ll - id) {
												jd = 0;
												id++;
											};
										};
										jc++;
										if (jc > lk - ic) {
											jc = 0;
											ic++;
										};
									};
									jb++;
									if (jb > lj - ib) {
										jb = 0;
										ib++;
									};
								};
								ja++;
								if (ja > li - ia) {
									ja = 0;
									ia++;
								};
							};

							/* switch to the Gamess holy order */
							if (gamess_order) {
								FullC<cdouble> full_trans_crt_gam(li, lj, lk, ll);
								full_trans_crt_gam.zero();
								GamessHolyOrder2E(full_trans_crt, full_trans_crt_gam, li, lj, lk, ll, true);
							};

							/* write down */
							full_trans_crt.load_pos_mx(posimx_crt_v[i], posimx_crt_v[j],
							                           posimx_crt_v[k], posimx_crt_v[l],
							                           i, j, k, l);
							#pragma omp critical
							{
								full_trans_crt.to_disk(tot_2E, keys.thrsh, ofs_2E);
							}

							///full_trans_crt.print();
						};
		}

		ofs_2E.close();
		cout << " Total number of two-electron integrals = " << tot_2E << endl;
	};

	/*quad fedr,fedi;
    bool flag;
    Faddeeva( 6.0, 3.0, fedr, fedi, flag );
    
    PrintQuadMath( fedr );
    PrintQuadMath( fedi );*/

	/* successful exit */
	t_end = omp_get_wtime();
	double duration = t_end - t_start;
	cout << endl;
	cout << " CPU time: " << setprecision(5) << fixed << duration << " s. " << endl;

	return EXIT_SUCCESS;
	/* */
};

void Init() {
	cout << endl;
	cout << "======================================================================================================" << endl;
	cout << "*                        PLANE WAVES + GAUSSIAN-TYPE ORBITALS INTEGRAL PACKAGE                       *" << endl;
	cout << "*                                 experimental version 0.2, Oct 2017                                 *" << endl;
	cout << "*                                  program written by Michal Lesiuk                                  *" << endl;
	cout << "*                                   at Quantum Chemistry Laboratory                                  *" << endl;
	cout << "*                          Faculty of Chemistry, University of Warsaw, Poland                        *" << endl;
	cout << "======================================================================================================" << endl;
	/* Thanks to J. Balcerzak, M. Jablczynska, and M. Szczygiel for their help at various stages of implementation and testing */
	cout << endl;
	cout.precision(15);
	cout << fixed;
	cout << scientific;

	belt[0] = 1;
	for (int i = 1; i < 500; i++)
		belt[i] = -belt[i - 1];

	fact[0] = 1.0;
	fact[1] = 1.0;
	for (int i = 2; i <= 160; i++)
		fact[i] = fact[i - 1] * i;

	for (int i = 0; i <= 100; i++)
		for (int j = 0; j <= i; j++) {
			binom[i][j] = fact[i] / (fact[j] * fact[i - j]);
			binom[j][i] = binom[i][j];
		};

	/*dfact[1] = 1.0;
    for(int i=0; i<=200; i+=2) dfact[i] = pow ( 2.0, i/2 ) * fact[i/2];
    for(int i=3; i<=200; i+=2) dfact[i] = dfact[i-2] * i;*/

	dfact[0] = 1.0;
	dfact[1] = 1.0;
	for (int i = 1; i < 200; i++)
		dfact[i + 1] = dfact[i] * (2 * i + 1);

	for (int i = 0; i <= 32; i++)
		for (int j = 0; j <= i; j++)
			omega[i][j] = sqrt((2.0 * i + 1) * fact[i - j] / (2.0 * fact[i + j]));

	int mm, pos, lx, ly, lz;
	for (int l = 0; l <= bas_lmax; l++)
		for (int m = -l; m <= l; m++) {
			mm = m + l;
			for (int i = 0; i <= l; i++)
				for (int j = 0; j <= l - i; j++) {
					pos = i * (2 * l + 3 - i) / 2 + j;
					lx = i;
					ly = j;
					lz = l - i - j;
					clmr[l][mm][pos] = CalcClmR(l, m, lx, ly, lz);
				};
		};

	int ia, ja, ka, shgi;
	for (int li = 0; li <= bas_lmax; li++) {
		shgi = crt_siz[li];

		ia = 0;
		ja = 0;
		for (int mi = 0; mi < shgi; mi++) {
			ka = li - ia - ja;
			xyz_norm[li][mi] = sqrt(dfact[ia] * dfact[ja] * dfact[ka]);

			ja++;
			if (ja > li - ia) {
				ja = 0;
				ia++;
			};
		};
	};
	/* */
};
