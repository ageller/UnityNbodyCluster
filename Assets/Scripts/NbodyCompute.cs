//https://gitlab.com/tzaeru/Nbody-simulation-tutorial-code-project
//https://www.tzaeru.com/site/compute-shader-basics-with-unity/
//23454.93
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class NbodyCompute : MonoBehaviour {

	public ComputeShader computeShader;
	public Transform cameraCenter;

	public int NumBodies = 30000;
	private int NumBodiesMax = 50000; //try setting up this system to start and then just limiting the max array location?

	public float velocityRange = 0.0f;
	public float Softening = 0.0f; //multiplied by the Earth particle radius
	//public float G = 6.67300e-11f;
	public float dvelMax = 0.0f; //max velocity to separate full vs. partial inelastic collisions
	public float rcolFac = 0.0f; //factor to multiply by the combined radii of particles to determine if a collision occurs
	public float restitution = 0.9f; //restitution coefficient for partially inelastic collisions
	public float rScale = 1.0f;//scaling for size of objects
	public float vScale = 1.0f;//scaling for velocity of bodies
	public float tScale = 1.0f;//scaling for dt (if tScale = 1, then time step is in Myr)
	public float neta = 0.5f;//neta for Reimer's mass loss
	public float Zmet = 0.02f;//metalicity (0.02 is solar)
	public float Rgc = 8000.0f; //parsecs, distance of cluster from center of galaxy
	public float hgc = 0.0f;//parsecs, height above the galactic plane
	public float ecc = 0.0f;//orbital eccentricity of galactic orbit (am I using this correctly?)
	public float DFractal = 3.0f;// fractal dimension =3.0 unfractal, =2.6 2/8 fractal, =2.0 4/8 fractal, =1.6 6/8 fractal, (2^D children per parent, following Goodwin & Whitworth 2004);
	//public float galMass = 1e11f; //solMass, mass of point-mass galaxy, mass inside Sun
	public float Rhm = 2.0f;//pc, half-mass radius (for cluster integration only)

	public static ComputeBuffer pos_buf;
	public static ComputeBuffer posCoM_buf;
	public static ComputeBuffer vel_buf;
	public static ComputeBuffer mass_buf;
	public static ComputeBuffer rad_buf;
	public static ComputeBuffer lum_buf;
	public static ComputeBuffer tgb_buf;
	public static ComputeBuffer lgb_buf;

	private float dt;
	private float time;
	private float fmGal;
	private float vesc;

	private Random.State state; //store the initial random seed

	public struct Cluster{
		public float[] pos_data;
		public float[] posCoM_data;
		public float[] vel_data;
		public float[] mass_data;
		public float[] rad_data;
		public float[] lum_data;
		public float[] tgb_data;
		public float[] lgb_data;
		public float Mcl;
	}

	public float[,] FractalStar;

	public Cluster cls;
	public Cluster cls0;

	public struct TL{
		public Vector4 times;
		public Vector4 lums;
	}

	private const int BLOCKSIZE = 128;
	//private const int floatSize = 12;//sizeof(float); = 4, but that doesn't seem correct! gives errors
	private const int floatSize = sizeof(float);// = 4, but that doesn't seem correct! gives errors

	private ButtonController bc;

	public GameObject galaxyCenter;

	public void tScaleReceiver(float val){
		tScale = Mathf.Pow(10.0f,val);
	}
	public void NumBodiesReceiver(int val){
		NumBodies = val;
		Random.state = state;
		doInit();
	}
	public void RhmReceiver(float val){
		Rhm = val;
		Random.state = state;
		doInit();
	}
	public void vScaleReceiver(float val){
		vScale = val;
		Random.state = state;
		doInit();
	}
	public void DFractalReceiver(float val){
		DFractal = val;
		Random.state = state;
		if (DFractal < 3.0f){
			Fractalize();
		}
		doInit();
	}
	public void RgcReceiver(float val){
		Rgc = val;
		Random.state = state;
		doInit();
	}
	public void eccReceiver(float val){
		ecc = val;
		Random.state = state;
		doInit();
	}
	public void doInit(){
		resetStars(NumBodies);
		clearStars(NumBodies);
		initPlummer(NumBodies);
		SetShaderData();

		InitCameraMover();
	}
	float KTGmass(){
		//Eq. 14 from Kroupa, Tout & Gilmore (1993)
		//Table 10, using alpha1 = 1.3 (not sure what is best)
		float c1 = 0.19f;
		float c2 = 1.55f;
		float c3 = 0.05f;
		float c4 = 0.6f;

		float X = Random.Range(0.0f, 1.0f);

		return 0.08f + (c1*Mathf.Pow(X, c2) + c3*Mathf.Pow(X, c4))/Mathf.Pow(1.0f - X, 0.58f);
	}


	float MWmass(float r){
		//Milky Way mass as a function of radius
		//NOTE: this assumes that the orbit is in the plane of the galactic disk.
		
		//Navarro-Frenk-White dark-matter halo mass as a function of radius for galaxy
		//https://en.wikipedia.org/wiki/Navarro%E2%80%93Frenk%E2%80%93White_profile
		//taking these values from table 1 in this paper: https://arxiv.org/pdf/1304.5127.pdf
		//https://ui.adsabs.harvard.edu/abs/2013JCAP...07..016N/abstract
		float Rs = 16100.0f;//pc Milky Way scale radius
		float rho0 = 0.014f ;//Msun/pc**3. Milky Way scale density
		//https://iopscience.iop.org/article/10.1086/589500/pdf
		//https://ui.adsabs.harvard.edu/abs/2008ApJ...684.1143X/abstract
		// float Rs = 22250.0f;//pc Milky Way scale radius
		// float rho0 = 0.0042f ;//Msun/pc**3. Milky Way scale density	
		float MNFW = 4.0f*CONSTANTS.pi*rho0*Mathf.Pow(Rs,3.0f)*(Mathf.Log((Rs + r)/Rs) - r/(Rs + r));

		//assume simple spherical ball of constant density for bulge
		float Mbulge = 1.5e10f;//Msun
		float Rbulge = 2000.0f;//pc
		float rhobulge = Mbulge/(4.0f/3.0f*CONSTANTS.pi*Mathf.Pow(Rbulge,3.0f));
		float MB = Mathf.Min(4.0f/3.0f*CONSTANTS.pi*Mathf.Pow(r, 3.0f)*rhobulge, Mbulge);

		//disk
		//also see here: https://academic.oup.com/mnras/article/366/3/899/993295
		float Rdisk = 3000.0f;//pc
		float Sdisk = 60.0f/Mathf.Exp(-8000/Rdisk);//Msun/pc**2
		float MD = 2.0f*CONSTANTS.pi*Sdisk*Rdisk*Rdisk*(1.0f - (1.0f + r/Rdisk)*Mathf.Exp(-r/Rdisk));

		return MNFW + MB + MD;
	}

	float MSRadFromMass(float m){
		//not sure where this came from
		float eta = 0.8f;

		if (m > 1.0f){ 
			eta = 0.57f;
		} 

		return Mathf.Pow(m, eta);

	}
	float MSLumFromMass(float m){
		//https://en.wikipedia.org/wiki/Mass%E2%80%93luminosity_relation
		float cons = 1.0f;
		float coeff = 4.0f;
		
		if (m < 0.43f){
			cons = 0.23f;
			coeff = 2.3f;
		}
		if (m > 2.0){// & m < 20.0){
			cons = 1.5f;
			coeff = 3.5f;
		}

		return cons*Mathf.Pow(m, coeff);
	}


	TL tlgb(float mass, float lMS){
		//times when the star reaches the base of the RGB, He burning, and AGB (Myr)
		//Hurley, Pols, Tout (2000), SSE paper
		//I would bet there are bugs here :)

		//Appendix A
		float psi = Mathf.Log(Zmet/0.02f);
		float a1 = Mathf.Pow(1.593890f, 3.0f) + Mathf.Pow(2.053038f, 3.0f)*psi + Mathf.Pow(1.231226f, 3.0f)*psi*psi + Mathf.Pow(2.327785f, 2.0f)*psi*psi*psi;
		float a2 = Mathf.Pow(2.706708f, 3.0f) + Mathf.Pow(1.483131f, 3.0f)*psi + Mathf.Pow(5.772723f, 2.0f)*psi*psi + Mathf.Pow(7.411230f, 1.0f)*psi*psi*psi;
		float a3 = Mathf.Pow(1.466143f, 2.0f) - Mathf.Pow(1.048442f, 2.0f)*psi - Mathf.Pow(6.795374f, 1.0f)*psi*psi - Mathf.Pow(1.391127f, 1.0f)*psi*psi*psi;
		float a4 = Mathf.Pow(4.141960f,-2.0f) + Mathf.Pow(4.564888f,-2.0f)*psi + Mathf.Pow(2.958542f,-2.0f)*psi*psi + Mathf.Pow(5.571483f,-3.0f)*psi*psi*psi;
		float a5 = Mathf.Pow(3.426349f,-1.0f);

		float b9  =  Mathf.Pow(2.751631f, 3.0f) + Mathf.Pow(3.557098f, 2.0f)*psi;
		float b10 = -Mathf.Pow(3.820831f,-2.0f) + Mathf.Pow(5.872664f,-2.0f)*psi;
		float b11p = Mathf.Pow(1.071738f, 2.0f) - Mathf.Pow(8.970339f, 1.0f)*psi - Mathf.Pow(3.949739f, 1.0f)*psi*psi;
		float b11 = b11p*b11p; 
		float b12 =  Mathf.Pow(7.348793f, 2.0f) - Mathf.Pow(1.531020f, 2.0f)*psi - Mathf.Pow(3.793700f, 2.0f)*psi*psi;
		float b13p = Mathf.Pow(9.219293f, 0.0f) - Mathf.Pow(2.005865f, 0.0f)*psi - Mathf.Pow(5.561309f,-1.0f)*psi*psi;
		float b13 = b13p*b13p;

		float b36 = Mathf.Pow(Mathf.Pow(1.445216f,-1.0f) - Mathf.Pow(6.180219f,-2.0f)*psi + Mathf.Pow(3.093878f,-2.0f)*psi*psi + Mathf.Pow(1.567090f,-2.0f)*psi*psi*psi, 4.0f);
		float b37 =     4.0f*(Mathf.Pow(1.304129f, 0.0f) + Mathf.Pow(1.395919f,-1.0f)*psi + Mathf.Pow(4.142455f,-3.0f)*psi*psi - Mathf.Pow(9.732503f,-3.0f)*psi*psi*psi);
		float b38 = Mathf.Pow(Mathf.Pow(5.114149f,-1.0f) - Mathf.Pow(1.160850f,-2.0f)*psi, 4.0f);

		float MHeF = 1.995f + 0.25f*psi + 0.087f*psi*psi; //eg. 2
		float MFGB = (13.048f*Mathf.Pow(Zmet/0.02f, 0.06f))/(1.0f + 0.0012f*Mathf.Pow(0.02f/Zmet, 1.27f)); //eq.3

		//////////////////////
		//RGB
		//time at the base of the rgb, eq. 4
		float tbrgb = (a1 + a2*Mathf.Pow(mass,4.0f) + a3*Mathf.Pow(mass, 5.5f) + Mathf.Pow(mass, 7.0f))/(a4*mass*mass + a5*Mathf.Pow(mass, 7.0f));
		//for now I am going to assume that the luminosity doesn't change on the MS, so that is an input to this function as the ZAMS lum
		//let's assume no evolution on the MS
		float lbrgb = lMS;

		//////////////////////
		//He burning
		//constants found on page 552 and 553 of the SSE paper, 
		float p = Mathf.Clamp(Mathf.Lerp(6.0f, 5.0f, (MHeF - mass)/(MHeF - 2.5f)), 5.0f, 6.0f);
		float D0 = 5.37f + 0.135f*psi;
		float D1 = Mathf.Max(Mathf.Max(-1.0f, 0.975f*D0 - 0.18f*mass), 0.5f*D0 - 0.06f*mass);
		float log10D = Mathf.Clamp(Mathf.Lerp(D0, D1, (MHeF - mass)/(MHeF - 2.5f)), D1, D0);
		float D = Mathf.Pow(10.0f, log10D);

		float AH = 1.44e-5f; //Msun / Lsun /Myr (this is the value on the RGB)
		//simplifying this and only taking the low-luminosity version
		float tinf1 = tbrgb + 1.0f/(AH*D*(p - 1.0f))*Mathf.Pow(D/lbrgb, (p - 1.0f)/p); //eq 40

		//eq. 49 for LHe
		float LHeI = (b11 + b12*Mathf.Pow(mass, 3.8f))/(b13 + mass*mass);
		if (mass < MHeF){
			float LHeMHeF = (b11 + b12*Mathf.Pow(MHeF, 3.8f))/(b13 + MHeF*MHeF);
			float alpha1 = b9*Mathf.Pow(MHeF, b10)/LHeMHeF - 1.0f;
			LHeI = b9*Mathf.Pow(mass, b10)/(1.0f + alpha1*Mathf.Exp(15.0f*(mass - MHeF))); 
		}

		// //eq. 51
		// float c = b17/Mathf.Pow(MFGB, 0.1f) + (b16*b17 - b14)/Mathf.Pow(MFGB, b15 + 0.1f);
		// float LminHe = LHeI*(b14 + c*Mathf.Pow(mass, b15+0.1))/(b16 + Mathf.Pow(mass, b15));

		// //eq. 53
		// float LZAHB = 
		//this looks OK, but what is "Lx"? only using low L equation.  Getting Lx requires a LOT more equations, and this THe bit is not working with out it (at least for the high mass stars, which makes sense).  I think this is not worth calculating! I will just do a check below to make sure that the time is between the tbrgb and tbagb
		float tHe = tinf1 - 1.0f/(AH*D*(p - 1.0f))*Mathf.Pow(D/LHeI, (p - 1.0f)/p); //eq 43


		//////////////////////
		//AGB 
		//I *think* this is what we want for the tabg, time that the star reaches the AGB.
		//this is when the TPABG (thermally pulsating AGB) begins
		AH = 7.66e-5f; //Msun / Lsun /Myr
		tinf1 = tbrgb + 1.0f/(AH*D*(p - 1.0f))*Mathf.Pow(D/lbrgb, (p - 1.0f)/p); //eq 40
		float Mcbagb = Mathf.Pow(b36*Mathf.Pow(mass, b37) + b38, 0.25f); //eq. 66
		float Lagb = D*Mathf.Pow(Mcbagb, p); //eq 31 (am I using this correctly?)
		float tbagb = tinf1 - 1.0f/(AH*D*(p - 1.0f))*Mathf.Pow(D/Lagb, (p - 1.0f)/p); //eq 70

		if (tHe > tbagb || tHe < tbrgb){
			//Debug.Log("BAD tHe : "+mass+" "+p+" "+D0+" "+D1+" "+log10D+" "+tinf1+" "+LHeI+" "+MHeF+" "+tHe);
			tHe = tbrgb + (tbagb - tbrgb)*0.1f;

		}
		
		//some time for the remnant phase
		float trem = tbagb + 1.0f;

		TL x;
		x.times = new Vector4(tbrgb, tHe, tbagb, trem);
		x.lums = new Vector4(lbrgb, LHeI, Lagb, 0.0f);

		return x;

	}



	//to test
	// float giantLum(float tbg, float lbg){
	// 	//giant luminosity vs. time
	// 	//eq. 25 in Jurley Pols & Tout 2000, the SSE paper

	// 	//simplifying this a little bit, in the paper these constants depend on mass and metallicity
	// 	float p = 5.0f;
	// 	float psi = 1.0f; //this could depend on metallicity
	// 	float D = Mathf.Pow(10.0f, 5.37f + 0.135f*psi);
	// 	float AH = 1.44e-5f; //Msun / Lsun /Myr

	// 	float tinf = tbg + 1.0f/(AH*D*(p - 1.0f))*Mathf.Pow(D/lbg, (p - 1.0f)/p);
	// 	float time = tbg;

	// 	return D*Mathf.Pow( (p - 1.0f)*AH*D*(tinf - time), p/(1.0f - p));
	// }
	// float giantRad(float m, float L){
	// 	//giant radius for solar mass stars (could use eq. 46 if we want metallicity dependence)
	// 	//eq. 48 in Hurley, Pols & Tout 2000

	// 	return 1.1f*Mathf.Pow(m, -0.3f)*(Mathf.Pow(L, 0.4f) + 0.383f*Mathf.Pow(L, 0.76f));
	// }


	// Use this for initialization or a circular orbit in the x-z plane
	//probably don't need Mcl, unless the galaxy mass is very small
	Vector3 CircularXY(float Mcl){
		
		//https://stackoverflow.com/questions/14845273/initial-velocity-vector-for-circular-orbit
		float xgc = Mathf.Sqrt(Rgc*Rgc + hgc*hgc);
		Vector3 dir = new Vector3(xgc, hgc, 0.0f);
		Vector3 dirNorm = dir/dir.magnitude;//(dir.x + dir.y + dir.z);
		float cosphi = dirNorm.x;
		float sinphi = dirNorm.z;
		float galMass = MWmass(Rgc);
		float initV = Mathf.Sqrt(CONSTANTS.G*(galMass + Mcl)/dir.magnitude);
		return new Vector3(-initV*sinphi, 0.0f, initV*cosphi); 


	}

	void InitStruct(int N){
		cls.pos_data = new float[N*3];
		cls.posCoM_data = new float[N*3];
		cls.vel_data = new float[N*3];
		cls.mass_data = new float[N*3];
		cls.rad_data = new float[N];
		cls.lum_data = new float[N];
		cls.tgb_data = new float[N*4];
		cls.lgb_data = new float[N*4];
		cls.Mcl = 0.0f;

		cls0.pos_data = new float[N*3];
		cls0.posCoM_data = new float[N*3];
		cls0.vel_data = new float[N*3];
		cls0.mass_data = new float[N*3];
		cls0.rad_data = new float[N];
		cls0.lum_data = new float[N];
		cls0.tgb_data = new float[N*4];
		cls0.lgb_data = new float[N*4];
		cls0.Mcl = 0.0f;

		FractalStar = new float[N, 7];

	}
	void clearStars(int N){
		//set these all to zero
		for (int i = N; i < NumBodiesMax; i++){

			cls.mass_data[i*3] = 0.0f; //solMass m(t)
			cls.mass_data[i*3 + 1] = 0.0f; //solMass m(t=0)
			cls.mass_data[i*3 + 2] = 0.0f; //solMass mcore
			cls.rad_data[i] = 0.0f; //solRad
			cls.lum_data[i] = 1.0f; //solar units

			cls.tgb_data[i*4] = 1e10f; //tbrgb Myr
			cls.tgb_data[i*4 + 1] = 1e10f; //tbHe Myr
			cls.tgb_data[i*4 + 2] = 1e10f; // tbagb Myr
			cls.tgb_data[i*4 + 3] = 1e10f; // trem Myr

			cls.lgb_data[i*4] = 1.0f; //tbrgb no evolution on MS for now
			cls.lgb_data[i*4 + 1] = 1.0f; //lbHe
			cls.lgb_data[i*4 + 2] = 1.0f; //lbAGB,
			cls.lgb_data[i*4 + 3] =  1.0f; //lrem (=0 for now)

			for (int j = 0; j < 7; j++) FractalStar[i,j] = 1.0f;

		}

	}

	void resetStars(int N){
		cls.Mcl = 0.0f;
		for (int i = 0; i < N; i++){
			cls.Mcl += cls0.mass_data[i*3];

			cls.mass_data[i*3] = cls0.mass_data[i*3]; //solMass m(t)
			cls.mass_data[i*3 + 1] = cls0.mass_data[i*3 + 1]; //solMass m(t=0)
			cls.mass_data[i*3 + 2] = cls0.mass_data[i*3 + 2]; //solMass mcore
			cls.rad_data[i] = cls0.rad_data[i]; //solRad
			cls.lum_data[i] = cls0.lum_data[i]; //solar units

			cls.tgb_data[i*4] = cls0.tgb_data[i*4]; //tbrgb Myr
			cls.tgb_data[i*4 + 1] = cls0.tgb_data[i*4 + 1]; //tbHe Myr
			cls.tgb_data[i*4 + 2] = cls0.tgb_data[i*4 + 2]; // tbagb Myr
			cls.tgb_data[i*4 + 3] = cls0.tgb_data[i*4 + 3]; // trem Myr

			cls.lgb_data[i*4] = cls0.lgb_data[i*4]; //tbrgb no evolution on MS for now
			cls.lgb_data[i*4 + 1] = cls0.lgb_data[i*4 + 1]; //lbHe
			cls.lgb_data[i*4 + 2] = cls0.lgb_data[i*4 + 2]; //lbAGB,
			cls.lgb_data[i*4 + 3] = cls0.lgb_data[i*4 + 3]; //lrem (=0 for now)

			//these may not be necessary to reset
			cls.pos_data[i*3] = cls0.pos_data[i*3];
			cls.pos_data[i*3 + 1] = cls0.pos_data[i*3 + 1];
			cls.pos_data[i*3 + 2] = cls0.pos_data[i*3 + 2];

			cls.posCoM_data[i*3] = cls0.posCoM_data[i*3];
			cls.posCoM_data[i*3 + 1] = cls0.posCoM_data[i*3 + 1];
			cls.posCoM_data[i*3 + 2] = cls0.posCoM_data[i*3 + 2];

			cls.vel_data[i*3] = cls0.vel_data[i*3];
			cls.vel_data[i*3 + 1] = cls0.vel_data[i*3 + 1];
			cls.vel_data[i*3 + 2] = cls0.vel_data[i*3 + 2];
		}
	}

	void initStars(int N){

		//first get all the masses from the IMF, and calculate the cluster mass
		//also set the radii and luminosities
		cls0.Mcl = 0.0f;
		for (int i = 0; i < N; i++){
			float m = KTGmass();
			float r = MSRadFromMass(m);
			float L = MSLumFromMass(m);

			cls0.Mcl += m;

			cls0.mass_data[i*3] = m; //solMass m(t)
			cls0.mass_data[i*3 + 1] = m; //solMass m(t=0)
			cls0.mass_data[i*3 + 2] = -1.0f; //solMass mcore
			cls0.rad_data[i] = r; //solRad
			cls0.lum_data[i] = L; //solar units

			TL t = tlgb(m, L);
			cls0.tgb_data[i*4] = Mathf.Clamp(t.times[0], 0.0f, 14000f); //tbrgb Myr
			cls0.tgb_data[i*4 + 1] = Mathf.Clamp(t.times[1], 0.0f, 14000f); //tbHe Myr
			cls0.tgb_data[i*4 + 2] = Mathf.Clamp(t.times[2], 0.0f, 14000f); // tbagb Myr
			cls0.tgb_data[i*4 + 3] = Mathf.Clamp(t.times[3], 0.0f, 14000f); // trem Myr

			cls0.lgb_data[i*4] = t.lums[0]; //tbrgb no evolution on MS for now
			cls0.lgb_data[i*4 + 1] = t.lums[1]; //lbHe
			cls0.lgb_data[i*4 + 2] = t.lums[2]; //lbAGB,
			cls0.lgb_data[i*4 + 3] = t.lums[3]; //lrem (=0 for now)

			//if (m > 5) Debug.Log(m+" "+L+" "+t.times[0]+" "+t.times[1]+" "+t.times[2]+" "+t.times[3]);
			//Debug.Log(m+" "+L+" "+cls.tgb_data[i*4]+" "+cls.tgb_data[i*4 + 1]+" "+cls.tgb_data[i*4 + 2]+" "+cls.tgb_data[i*4 + 3]);
		}


	}

	//converted from McCluster by Andreas Kupper
	//forcing spherical symmetry, and setting to "radial"
	float drand48(){
		//to simplify the conversion
		return Random.Range(0.0f, 1.0f);
	}
	float get_gauss(){
		float[] random = new float[2];
		float p, q;
		q = 2.0f;
		while (q > 1.0f){
			random[0] = 2.0f*drand48()-1.0f; 
			random[1] = 2.0f*drand48()-1.0f; 
			q = random[0]*random[0]+random[1]*random[1];
		}
		
		p = Mathf.Sqrt(-2.0f*Mathf.Log(q)/q);
		return random[0]*p;
		
	}

	void Fractalize(){
	
		int N = NumBodies;

		bool radial = true;
		bool symmetry = true;

		int i = 0;
		int j = 0;
		int h = 0;
		int Nparent = 0;
		int Nparentlow = 0;
		int Ntot = (int)Mathf.Round(128.0f*Mathf.Pow(8f,Mathf.Ceil(Mathf.Log(NumBodies)/Mathf.Log(8f))));
		int Ntotorg = Ntot;
		float l = 2.0f;
		float prob = Mathf.Pow(2.0f, DFractal - 3.0f);
		Debug.Log("D "+DFractal+" "+prob);
		float scatter = 0.1f;
		if (radial) scatter = 0.01f; 
		float vx = 0.0f;
		float vy = 0.0f;
		float vz = 0.0f;
		float vscale = 0.0f;
		int subi = 0;
		float morescatter = 0.1f; //0.1 looks good


		//temporary array for fractalized structure
		float[,] star_temp = new float[Ntot, 7]; //not sure why he's using 7 when there are only 6 positions needed?

		Nparent = 0; //ur-star
		Nparentlow = Nparent;
		star_temp[Nparent,1] = 0.0f;//x
		star_temp[Nparent,2] = 0.0f;//y
		star_temp[Nparent,3] = 0.0f;//z
		star_temp[Nparent,4] = 0.0f;//vx
		star_temp[Nparent,5] = 0.0f;//vy
		star_temp[Nparent,6] = 0.0f;//vz
		Nparent++;
		
		
		while (Nparent+i*8<Ntot) {
			l /= 2.0f;
			i = 0;
			for (;Nparentlow<Nparent;Nparentlow++) {
				subi = 0;

				if ((drand48()<prob && Nparent+i<Ntot) || ((Nparent == 1) && (symmetry))) {
					star_temp[Nparent+i,0] = 1.0f;
					star_temp[Nparent+i,1] = star_temp[Nparentlow,1] + l/2.0f + l*scatter*get_gauss();
					star_temp[Nparent+i,2] = star_temp[Nparentlow,2] + l/2.0f + l*scatter*get_gauss();
					star_temp[Nparent+i,3] = star_temp[Nparentlow,3] + l/2.0f + l*scatter*get_gauss();
					star_temp[Nparent+i,4] = get_gauss();
					star_temp[Nparent+i,5] = get_gauss();
					star_temp[Nparent+i,6] = get_gauss();
					i++;
					subi++;
				}
				if ((drand48()<prob && Nparent+i<Ntot) || ((Nparent == 1) && (symmetry))) {
					star_temp[Nparent+i,0] = 1.0f;
					star_temp[Nparent+i,1] = star_temp[Nparentlow,1] + l/2.0f + l*scatter*get_gauss();
					star_temp[Nparent+i,2] = star_temp[Nparentlow,2] + l/2.0f + l*scatter*get_gauss();
					star_temp[Nparent+i,3] = star_temp[Nparentlow,3] - l/2.0f + l*scatter*get_gauss();
					star_temp[Nparent+i,4] = get_gauss();
					star_temp[Nparent+i,5] = get_gauss();
					star_temp[Nparent+i,6] = get_gauss();
					i++;
					subi++;
				}
				if ((drand48()<prob && Nparent+i<Ntot) || ((Nparent == 1) && (symmetry))) {
					star_temp[Nparent+i,0] = 1.0f;
					star_temp[Nparent+i,1] = star_temp[Nparentlow,1] + l/2.0f + l*scatter*get_gauss();
					star_temp[Nparent+i,2] = star_temp[Nparentlow,2] - l/2.0f + l*scatter*get_gauss();
					star_temp[Nparent+i,3] = star_temp[Nparentlow,3] + l/2.0f + l*scatter*get_gauss();
					star_temp[Nparent+i,4] = get_gauss();
					star_temp[Nparent+i,5] = get_gauss();
					star_temp[Nparent+i,6] = get_gauss();
					i++;
					subi++;
				}
				if ((drand48()<prob && Nparent+i<Ntot) || ((Nparent == 1) && (symmetry))) {
					star_temp[Nparent+i,0] = 1.0f;
					star_temp[Nparent+i,1] = star_temp[Nparentlow,1] - l/2.0f + l*scatter*get_gauss();
					star_temp[Nparent+i,2] = star_temp[Nparentlow,2] + l/2.0f + l*scatter*get_gauss();
					star_temp[Nparent+i,3] = star_temp[Nparentlow,3] + l/2.0f + l*scatter*get_gauss();
					star_temp[Nparent+i,4] = get_gauss();
					star_temp[Nparent+i,5] = get_gauss();
					star_temp[Nparent+i,6] = get_gauss();
					i++;
					subi++;
				}
				if ((drand48()<prob && Nparent+i<Ntot) || ((Nparent == 1) && (symmetry))) {
					star_temp[Nparent+i,0] = 1.0f;
					star_temp[Nparent+i,1] = star_temp[Nparentlow,1] + l/2.0f + l*scatter*get_gauss();
					star_temp[Nparent+i,2] = star_temp[Nparentlow,2] - l/2.0f + l*scatter*get_gauss();
					star_temp[Nparent+i,3] = star_temp[Nparentlow,3] - l/2.0f + l*scatter*get_gauss();
					star_temp[Nparent+i,4] = get_gauss();
					star_temp[Nparent+i,5] = get_gauss();
					star_temp[Nparent+i,6] = get_gauss();
					i++;
					subi++;
				}
				if ((drand48()<prob && Nparent+i<Ntot) || ((Nparent == 1) && (symmetry))) {
					star_temp[Nparent+i,0] = 1.0f;
					star_temp[Nparent+i,1] = star_temp[Nparentlow,1] - l/2.0f + l*scatter*get_gauss();
					star_temp[Nparent+i,2] = star_temp[Nparentlow,2] + l/2.0f + l*scatter*get_gauss();
					star_temp[Nparent+i,3] = star_temp[Nparentlow,3] - l/2.0f + l*scatter*get_gauss();
					star_temp[Nparent+i,4] = get_gauss();
					star_temp[Nparent+i,5] = get_gauss();
					star_temp[Nparent+i,6] = get_gauss();
					i++;
					subi++;
				}
				if ((drand48()<prob && Nparent+i<Ntot) || ((Nparent == 1) && (symmetry))) {
					star_temp[Nparent+i,0] = 1.0f;
					star_temp[Nparent+i,1] = star_temp[Nparentlow,1] - l/2.0f + l*scatter*get_gauss();
					star_temp[Nparent+i,2] = star_temp[Nparentlow,2] - l/2.0f + l*scatter*get_gauss();
					star_temp[Nparent+i,3] = star_temp[Nparentlow,3] + l/2.0f + l*scatter*get_gauss();
					star_temp[Nparent+i,4] = get_gauss();
					star_temp[Nparent+i,5] = get_gauss();
					star_temp[Nparent+i,6] = get_gauss();
					i++;
					subi++;
				}
				if ((drand48()<prob && Nparent+i<Ntot) || ((Nparent == 1) && (symmetry))) {
					star_temp[Nparent+i,0] = 1.0f;
					star_temp[Nparent+i,1] = star_temp[Nparentlow,1] - l/2.0f + l*scatter*get_gauss();
					star_temp[Nparent+i,2] = star_temp[Nparentlow,2] - l/2.0f + l*scatter*get_gauss();
					star_temp[Nparent+i,3] = star_temp[Nparentlow,3] - l/2.0f + l*scatter*get_gauss();
					star_temp[Nparent+i,4] = get_gauss();
					star_temp[Nparent+i,5] = get_gauss();
					star_temp[Nparent+i,6] = get_gauss();
					i++;
					subi++;
				}
				
				//re-scaling of sub-group
				if (subi != 0) {
					vx = 0.0f; vy = 0.0f; vz = 0.0f;
					vscale = 0.0f;
					for (j=Nparent+i-subi;j<Nparent+i;j++) {
						vx += star_temp[j,4];
						vy += star_temp[j,5];
						vz += star_temp[j,6];
					}
					vx /= 1.0f*subi;
					vy /= 1.0f*subi;
					vz /= 1.0f*subi;
					for (j=Nparent+i-subi;j<Nparent+i;j++) {
						star_temp[j,4] -= vx;
						star_temp[j,5] -= vy;
						star_temp[j,6] -= vz;					
					}
					if (subi-1 != 0) {
						for (j=Nparent+i-subi;j<Nparent+i;j++) {
							vscale += star_temp[j,4]*star_temp[j,4]+star_temp[j,5]*star_temp[j,5]+star_temp[j,6]*star_temp[j,6];
						}
						vscale = Mathf.Sqrt(vscale/(subi-1));
					} else {
						vscale = 1.0f;
					}
					for (j=Nparent+i-subi;j<Nparent+i;j++) {
						star_temp[j,4] = star_temp[j,4]/vscale + star_temp[Nparentlow,4]; //add bulk velocity of parent
						star_temp[j,5] = star_temp[j,5]/vscale + star_temp[Nparentlow,5];
						star_temp[j,6] = star_temp[j,6]/vscale + star_temp[Nparentlow,6];
					}				
				}
				
			}
			Nparent+=i;
		}
		Ntot = Nparent;

		float[] cmr = new float[7];//centre-of-mass correction
		for (j=0; j<7; j++) cmr[j] = 0.0f;
		
		for (j=0; j<Ntot; j++) {
			for (i=1;i<7;i++) 
				cmr[i] += star_temp[j,i]; 
		} 
		
		for (j=0; j<Ntot; j++) {
			for (i=1;i<7;i++)
				star_temp[j,i] -= cmr[i]/Ntot;
		}
			
		i = 0;//randomly select stars from sample with r < 1.0
		for (i=0; i<N; i++) {
			do{ 
				j = (int)Mathf.Round(drand48()*Ntot);
			} while (star_temp[j,0] == 0.0f || (Mathf.Sqrt(Mathf.Pow(star_temp[j,1],2)+Mathf.Pow(star_temp[j,2],2)+Mathf.Pow(star_temp[j,3],2)) > 1.0)); 
			for (h=1;h<7;h++) FractalStar[i,h] = star_temp[j,h];
			star_temp[j,0] = 0.0f;
		}

		float r, r_norm;
		float vnorm = 0.0f;
		float[] start = new float[4];
		if (radial) {
			for (i=0;i<N;i++) {
				vnorm += Mathf.Sqrt(Mathf.Pow(FractalStar[i,4],2)+Mathf.Pow(FractalStar[i,5],2)+Mathf.Pow(FractalStar[i,6],2));
			}
			vnorm /= N;
			//this for loop was not in McCluster... but it raises an error here because i = N (and the star array is only of size N)
			//I added the loop, which I think makes sense?
			for (i=0;i<N;i++) {
				for (h=4;h<7;h++) FractalStar[i,h] /= 0.5f*vnorm;
			}
			
			for (i=0;i<N;i++) {
				r = Mathf.Sqrt(Mathf.Pow(FractalStar[i,1],2)+Mathf.Pow(FractalStar[i,2],2)+Mathf.Pow(FractalStar[i,3],2));
				r_norm = Mathf.Pow(r,3);
				do{ 
					for (h=1;h<4;h++) start[h] = FractalStar[i,h]*r_norm/r + Mathf.Pow(r_norm/r,3)*morescatter*get_gauss();
				} while (Mathf.Sqrt(Mathf.Pow(start[1],2)+Mathf.Pow(start[2],2)+Mathf.Pow(start[3],2))>1.0);
				for (h=1;h<4;h++) FractalStar[i,h] = start[h];
			}
		}
		//for (j=0;j<Ntotorg;j++) free (star_temp[j]);
		//free(star_temp);
				
	}

	void initUnscaledPlummer(int N){
		//For rendering, I need this to be in the cluster's center-of-mass frame, but how can I do that??
		//then in the compute shader I need to have the galaxy orbiting 

		int i;

		float rpl = 1.0f;

		//now set the positions and velocities
		for (i = 0; i < N; i++){

			//positions
			float X1 = Random.Range(0.0f, 1.0f);
			float psi = Mathf.Pow(Mathf.Pow(X1, -2.0f/3.0f) - 1.0f, -0.5f);
			float r = psi*rpl;

			float X2 = Random.Range(0.0f, 1.0f);
			float X3 = Random.Range(0.0f, 2.0f*CONSTANTS.pi);
			float z = (2.0f*X2 - 1.0f)*r;
			float xy = Mathf.Sqrt(r*r - z*z);
			float x = xy*Mathf.Cos(X3);
			float y = xy*Mathf.Sin(X3);

			cls0.pos_data[i*3] = x ;
			cls0.pos_data[i*3 + 1] = y;
			cls0.pos_data[i*3 + 2] = z;

			cls0.posCoM_data[i*3] = x;
			cls0.posCoM_data[i*3 + 1] = y;
			cls0.posCoM_data[i*3 + 2] = z;

			//not sure I can give unscaled velocities?
			cls0.vel_data[i*3] = 0.0f;
			cls0.vel_data[i*3 + 1] = 0.0f;
			cls0.vel_data[i*3 + 2] = 0.0f;
			
		}

	}
	
	void initPlummer(int N){
		//this assumes that an unscaled plummer position array is available in cls0

		int i;

		float xgc = Mathf.Sqrt(Rgc*Rgc + hgc*hgc);
		Vector3 center = new Vector3(xgc, hgc, 0.0f);
		float rpl = Rhm*Mathf.Sqrt(Mathf.Pow(2.0f, 2.0f/3.0f) -1.0f);
		Vector3 vGc = CircularXY(cls.Mcl)*(1.0f - ecc); //is this correct?
		float galMass = MWmass(Rgc);
		//from McCluster
		float Rtide = Mathf.Pow(1.0f*cls.Mcl/(3.0f*galMass), 1.0f/3.0f)*Rgc;

		Debug.Log("cluster mass, velocity, Rhm, rtide, vscale, Mgal = "+cls.Mcl+" "+vGc+" "+Rhm+" "+Rtide+" "+vScale+" "+galMass);

		//now set the positions and velocities
		for (i = 0; i < N; i++){

			//positions

			float x = cls0.pos_data[i*3]*rpl;
			float y = cls0.pos_data[i*3 + 1]*rpl;
			float z = cls0.pos_data[i*3 + 2]*rpl;
			if (DFractal < 3){
				//apply the fractal pattern?, not sure how to scale this properly
				x = FractalStar[i,1]*rpl;
				y = FractalStar[i,2]*rpl;
				z = FractalStar[i,3]*rpl;
			}
			float r = Mathf.Sqrt(x*x + y*y + z*z);
			float psi = r/rpl;

			//limit by approximate tidal radius (this will place them in the Plummer profile even if Fractal is desired...)
			if (r > Rtide){
				while (r > Rtide){
					float X1 = Random.Range(0.0f, 1.0f);
					psi = Mathf.Pow(Mathf.Pow(X1, -2.0f/3.0f) - 1.0f, -0.5f);
					r = psi*rpl;
				}

				float X2 = Random.Range(0.0f, 1.0f);
				float X3 = Random.Range(0.0f, 2.0f*CONSTANTS.pi);
				z = (2.0f*X2 - 1.0f)*r;
				float xy = Mathf.Sqrt(r*r - z*z);
				x = xy*Mathf.Cos(X3);
				y = xy*Mathf.Sin(X3);
			}




			cls.pos_data[i*3] = x + center.x;
			cls.pos_data[i*3 + 1] = y + center.y;
			cls.pos_data[i*3 + 2] = z + center.z;

			cls.posCoM_data[i*3] = x;
			cls.posCoM_data[i*3 + 1] = y;
			cls.posCoM_data[i*3 + 2] = z;

			//velocities
			float X4 = Random.Range(0.0f, 1.0f);
			float X5 = Random.Range(0.0f, 2.0f*CONSTANTS.pi);
			//velocities
			//It seems like there is a factor of 2 missing somewhere.  Using these velocities, the cluster always flies apart
			float phi = -CONSTANTS.G*cls.Mcl/rpl*Mathf.Pow(1.0f + Mathf.Pow(r/rpl, 2.0f), -0.5f);
			vesc = Mathf.Sqrt(-2.0f*phi);
			float v = psi*vesc*vScale*0.5f; //added factor of 0.5 here, not sure it's correct, but it seems needed to avoid the cluster exploding initially
			float vz = (2.0f*X4 - 1.0f)*v;
			float vxy = Mathf.Sqrt(v*v - vz*vz);
			float vx = vxy*Mathf.Cos(X5);
			float vy = vxy*Mathf.Sin(X5);

			//I don't think there's anything special about the velocities in Fractalize (rigtht?)
			// if (DFractal < 3){
			// 	//apply the fractal pattern?, not sure how to scale this properly
			// 	vx = FractalStar[i,4];
			// 	vy = FractalStar[i,5];
			// 	vz = FractalStar[i,6];
			// }
			cls.vel_data[i*3] = vx + vGc.x;
			cls.vel_data[i*3 + 1] = vy + vGc.y;
			cls.vel_data[i*3 + 2] = vz + vGc.z;
			
		}



		//try this for the center of mass calculation?
		i = N - 1;
		cls.pos_data[i*3] = center.x;
		cls.pos_data[i*3 + 1] = center.y;
		cls.pos_data[i*3 + 2] = center.z;
		cls.vel_data[i*3] = vGc.x;
		cls.vel_data[i*3 + 1] = vGc.y;
		cls.vel_data[i*3 + 2] = vGc.z;
		cls.mass_data[i*3] = cls.Mcl; 
		cls.rad_data[i] = 0.0f; 

	}

	void initPlummerOrg(int N){
		//doesn't require a precomputed Plummer position array

		int i;

		float xgc = Mathf.Sqrt(Rgc*Rgc + hgc*hgc);
		Vector3 center = new Vector3(xgc, hgc, 0.0f);
		float rpl = Rhm*Mathf.Sqrt(Mathf.Pow(2.0f, 2.0f/3.0f) -1.0f);
		Vector3 vcirc = CircularXY(cls.Mcl);
		float galMass = MWmass(Rgc);
		//from McCluster
		float Rtide = Mathf.Pow(1.0f*cls.Mcl/(3.0f*galMass), 1.0f/3.0f)*Rgc;
		fmGal = CONSTANTS.G*galMass/(Rgc*Rgc);

		Debug.Log("cluster mass, velocity, Rhm, rtide, vscale = "+cls.Mcl+" "+vcirc+" "+Rhm+" "+Rtide+" "+vScale);

		//now set the positions and velocities
		for (i = 0; i < N; i++){

			//positions
			float X1 = Random.Range(0.0f, 1.0f);
			float psi = Mathf.Pow(Mathf.Pow(X1, -2.0f/3.0f) - 1.0f, -0.5f);
			float r = psi*rpl;
			//limit by approximate tidal radius, but this breaks the re-randomizes the positions
			// while (r > Rtide){
			// 	X1 = Random.Range(0.0f, 1.0f);
			// 	psi = Mathf.Pow(Mathf.Pow(X1, -2.0f/3.0f) - 1.0f, -0.5f);
			// 	r = psi*rpl;
			// }

			float X2 = Random.Range(0.0f, 1.0f);
			float X3 = Random.Range(0.0f, 2.0f*CONSTANTS.pi);
			float z = (2.0f*X2 - 1.0f)*r;
			float xy = Mathf.Sqrt(r*r - z*z);
			float x = xy*Mathf.Cos(X3);
			float y = xy*Mathf.Sin(X3);

			float X4 = Random.Range(0.0f, 1.0f);
			float X5 = Random.Range(0.0f, 2.0f*CONSTANTS.pi);
			//velocities
			//It seems like there is a factor of 2 missing somewhere.  Using these velocities, the cluster always flies apart
			float phi = -CONSTANTS.G*cls.Mcl/rpl*Mathf.Pow(1.0f + Mathf.Pow(r/rpl, 2.0f), -0.5f);
			vesc = Mathf.Sqrt(-2.0f*phi);
			float v = psi*vesc*vScale*0.5f; //added factor of 0.5 here, not sure it's correct, but it seems needed to avoid collapse
			float vz = (2.0f*X4 - 1.0f)*v;
			float vxy = Mathf.Sqrt(v*v - vz*vz);
			float vx = vxy*Mathf.Cos(X5);
			float vy = vxy*Mathf.Sin(X5);

			cls.pos_data[i*3] = x + center.x;
			cls.pos_data[i*3 + 1] = y + center.y;
			cls.pos_data[i*3 + 2] = z + center.z;

			cls.posCoM_data[i*3] = x;
			cls.posCoM_data[i*3 + 1] = y;
			cls.posCoM_data[i*3 + 2] = z;

			cls.vel_data[i*3] = vx + vcirc.x;
			cls.vel_data[i*3 + 1] = vy + vcirc.y;
			cls.vel_data[i*3 + 2] = vz + vcirc.z;
			
		}



		//try this for the center of mass calculation?
		i = N - 1;
		cls.pos_data[i*3] = center.x;
		cls.pos_data[i*3 + 1] = center.y;
		cls.pos_data[i*3 + 2] = center.z;
		cls.vel_data[i*3] = vcirc.x;
		cls.vel_data[i*3 + 1] = vcirc.y;
		cls.vel_data[i*3 + 2] = vcirc.z;
		cls.mass_data[i*3] = cls.Mcl; 
		cls.rad_data[i] = 0.0f; 

	}
	void InitBuffers(int N){
		//Random.seed = seed;

		//initialize all the buffers
		//something is strange about the amount of size needed here;
		// with floatSize = 12, this works for 10k, 15k and 30k particles, but not 20k
		pos_buf = new ComputeBuffer(N, 3*floatSize, ComputeBufferType.Default);
		posCoM_buf = new ComputeBuffer(N, 3*floatSize, ComputeBufferType.Default);
		vel_buf = new ComputeBuffer(N, 3*floatSize, ComputeBufferType.Default);
		mass_buf = new ComputeBuffer(N, 3*floatSize, ComputeBufferType.Default); //mnow, m0, mcore (currently don't have mcore)
		rad_buf = new ComputeBuffer(N, floatSize, ComputeBufferType.Default);
		lum_buf = new ComputeBuffer(N, floatSize, ComputeBufferType.Default);
		tgb_buf = new ComputeBuffer(N, 4*floatSize, ComputeBufferType.Default);
		lgb_buf = new ComputeBuffer(N, 4*floatSize, ComputeBufferType.Default);

		// These global buffers apply to every shader with these buffers defined.
		Shader.SetGlobalBuffer(Shader.PropertyToID("position"), pos_buf);
		Shader.SetGlobalBuffer(Shader.PropertyToID("positionCoM"), posCoM_buf);
		Shader.SetGlobalBuffer(Shader.PropertyToID("velocity"), vel_buf);
		Shader.SetGlobalBuffer(Shader.PropertyToID("mass"), mass_buf);
		Shader.SetGlobalBuffer(Shader.PropertyToID("radius"), rad_buf);
		Shader.SetGlobalBuffer(Shader.PropertyToID("luminosity"), lum_buf);
		Shader.SetGlobalBuffer(Shader.PropertyToID("tgb"), tgb_buf);
		Shader.SetGlobalBuffer(Shader.PropertyToID("lgb"), lgb_buf);

		//For rendering, I need this to be in the cluster's center-of-mass frame, but how can I do that??
		//then in the compute shader I need to have the galaxy orbiting 

	}

	void SetShaderData(){
		pos_buf.SetData(cls.pos_data);
		posCoM_buf.SetData(cls.posCoM_data);
		vel_buf.SetData(cls.vel_data);
		mass_buf.SetData(cls.mass_data);
		rad_buf.SetData(cls.rad_data);
		lum_buf.SetData(cls.lum_data);
		tgb_buf.SetData(cls.tgb_data);
		lgb_buf.SetData(cls.lgb_data);
	}

	void DisposeOfBuffers(){
		pos_buf.Dispose();
		posCoM_buf.Dispose();
		vel_buf.Dispose();
		mass_buf.Dispose();
		rad_buf.Dispose();
		lum_buf.Dispose();
		tgb_buf.Dispose();
		lgb_buf.Dispose();
	}


	void Start(){
		bc = GameObject.Find("ButtonController").GetComponent<ButtonController>();

		InitBuffers(NumBodiesMax);
		InitStruct(NumBodiesMax);

		initStars(NumBodiesMax); //masses, radii, and other SEV items
		initUnscaledPlummer(NumBodiesMax);

		//apply the fractal structure?
		if (DFractal < 3.0f){
			Fractalize();
		}

		//to save the initial state
		state = Random.state; //for some reason this still looks slightly different than after using the sliders once??

		//this will copy from the initial state, and then (re)do the Plummer positions and velocities, and set the Shader data
		doInit();
	}
	


	// FixedUpdate is called many times per frame, Update is called once per frame
	// void FixedUpdate(){
	// 	if (!bc.paused)	Time.timeScale = tScale;
	// 	dt = Time.fixedDeltaTime;
	void Update(){

		//sending units of pc, MSun, Myr, but sizes in RSun
		if (!bc.paused)	Time.timeScale = tScale;
		dt = Time.deltaTime;
		time += dt;//*tScale;
		computeShader.SetFloat("deltaTime", dt);
		computeShader.SetFloat("time", time);
		computeShader.SetInt("NumBodies", NumBodies);
		computeShader.SetFloat("Softening", Softening);
		computeShader.SetFloat("G", CONSTANTS.G);
		computeShader.SetFloat("dvelMax",dvelMax); 
		computeShader.SetFloat("rcolFac",rcolFac); 
		computeShader.SetFloat("restitution", restitution);
		computeShader.SetFloat("neta", neta);

		//Debug.Log(dt+" "+Time.timeScale+" "+bc.paused);

		var numberOfGroups = Mathf.CeilToInt((float) NumBodies/BLOCKSIZE);
		computeShader.Dispatch(computeShader.FindKernel("CoMCompute"), 1, 1, 1);
		computeShader.Dispatch(computeShader.FindKernel("CSMain"), numberOfGroups, 1, 1);
	}

	void OnDestroy(){
		DisposeOfBuffers();

	}


	//////////////////////////////////
	//to move the camera so that the cluster actually orbits the galaxy when rendered
	void Attract()
	{

		float epsilon = 0.0f;//softening

		//to otherBody
		Vector3 direction = cameraCenter.GetComponent<Rigidbody>().position; //galaxy center is fixed at 0,0,0
		float distance = direction.magnitude;
		float galMass = MWmass(distance);

		float forceMagnitude = -CONSTANTS.G*(cls.Mcl*galMass)/(Mathf.Pow(distance, 2) + Mathf.Pow(epsilon,2));
		Vector3 force = direction.normalized*forceMagnitude;

		cameraCenter.GetComponent<Rigidbody>().AddForce(force);
	}

	void InitCameraMover(){

		float xgc = Mathf.Sqrt(Rgc*Rgc + hgc*hgc);
		Vector3 pos = new Vector3(xgc, hgc, 0.0f);

		cameraCenter.GetComponent<Transform>().position = pos;
		cameraCenter.GetComponent<Rigidbody>().position = pos;
		cameraCenter.GetComponent<Rigidbody>().mass = cls.Mcl;
		cameraCenter.GetComponent<Rigidbody>().velocity = CircularXY(cls.Mcl)*(1.0f - ecc); //am I using this correctly?

	}

	void FixedUpdate(){
	//void Update(){
		Attract();
	}
}
