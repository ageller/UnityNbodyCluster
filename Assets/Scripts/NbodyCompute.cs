//https://gitlab.com/tzaeru/Nbody-simulation-tutorial-code-project
//https://www.tzaeru.com/site/compute-shader-basics-with-unity/
//23454.93
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class NbodyCompute : MonoBehaviour {

	public ComputeShader shader;

	public int NumBodies = 30000;

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
	public float Rgc = 8500.0f; //parsecs, distance of cluster from center of galaxy
	public float galMass = 1e11f; //solMass, mass of point-mass galaxy, mass inside Sun

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

	public struct Particle{
		public float[] pos_data;
		public float[] posCoM_data;
		public float[] vel_data;
		public float[] mass_data;
		public float[] rad_data;
		public float[] lum_data;
		public float[] tgb_data;
		public float[] lgb_data;
	}


	public struct TL{
		public Vector4 times;
		public Vector4 lums;
	}

	private const int BLOCKSIZE = 128;
	//private const int floatSize = 12;//sizeof(float); = 4, but that doesn't seem correct! gives errors
	private const int floatSize = sizeof(float);// = 4, but that doesn't seem correct! gives errors

	private CameraMover cameraMover;
	private ButtonController bc;

	public GameObject galaxyCenter;

	public void tScaleReceiver(float val){
		tScale = Mathf.Pow(10.0f,val);
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
		Vector3 dir = new Vector3(Rgc, 0.0f, 0.0f);
		Vector3 dirNorm = dir/dir.magnitude;//(dir.x + dir.y + dir.z);
		float cosphi = dirNorm.x;
		float sinphi = dirNorm.z;

		float initV = Mathf.Sqrt(CONSTANTS.G*(galMass + Mcl)/dir.magnitude);
		return new Vector3(-initV*sinphi, 0.0f, initV*cosphi); 


	}

	Particle initPlummer(Vector3 center, float Rhm, int N0, int Nf){
		//For rendering, I need this to be in the cluster's center-of-mass frame, but how can I do that??
		//then in the compute shader I need to have the galaxy orbiting 

		int i;

		Particle p;
		p.pos_data = new float[NumBodies*3];
		p.posCoM_data = new float[NumBodies*3];
		p.vel_data = new float[NumBodies*3];
		p.mass_data = new float[NumBodies*3];
		p.rad_data = new float[NumBodies];
		p.lum_data = new float[NumBodies];
		p.tgb_data = new float[NumBodies*4];
		p.lgb_data = new float[NumBodies*4];


		float rpl = Rhm*Mathf.Sqrt(Mathf.Pow(2.0f, 2.0f/3.0f) -1.0f);

		//first get all the masses from the IMF, and calculate the cluster mass
		//also set the radii and luminosities
		float Mcl = 0.0f;
		for (i = N0; i < Nf; i++){
			float m = KTGmass();
			float r = MSRadFromMass(m);
			float L = MSLumFromMass(m);

			Mcl += m;

			p.mass_data[i*3] = m; //solMass m(t)
			p.mass_data[i*3 + 1] = m; //solMass m(t=0)
			p.mass_data[i*3 + 2] = -1.0f; //solMass mcore
			p.rad_data[i] = r; //solRad
			p.lum_data[i] = L; //solar units

			TL t = tlgb(m, L);
			p.tgb_data[i*4] = Mathf.Clamp(t.times[0], 0.0f, 14000f); //tbrgb Myr
			p.tgb_data[i*4 + 1] = Mathf.Clamp(t.times[1], 0.0f, 14000f); //tbHe Myr
			p.tgb_data[i*4 + 2] = Mathf.Clamp(t.times[2], 0.0f, 14000f); // tbagb Myr
			p.tgb_data[i*4 + 3] = Mathf.Clamp(t.times[3], 0.0f, 14000f); // trem Myr

			p.lgb_data[i*4] = t.lums[0]; //tbrgb no evolution on MS for now
			p.lgb_data[i*4 + 1] = t.lums[1]; //lbHe
			p.lgb_data[i*4 + 2] = t.lums[2]; //lbAGB,
			p.lgb_data[i*4 + 3] = t.lums[3]; //lrem (=0 for now)


			//if (m > 5) Debug.Log(m+" "+L+" "+t.times[0]+" "+t.times[1]+" "+t.times[2]+" "+t.times[3]);
			//Debug.Log(m+" "+L+" "+p.tgb_data[i*4]+" "+p.tgb_data[i*4 + 1]+" "+p.tgb_data[i*4 + 2]+" "+p.tgb_data[i*4 + 3]);
		}
		Vector3 vcirc = CircularXY(Mcl);

		Debug.Log("cluster mass = "+Mcl);
		Debug.Log("cluster velocity = "+vcirc);
		cameraMover.SendMessage("MclReceiver", Mcl);
		fmGal = CONSTANTS.G*galMass/(Rgc*Rgc);

		//now set the positions and velocities
		for (i = N0; i < Nf; i++){
			float X1 = Random.Range(0.0f, 1.0f);
			float X2 = Random.Range(0.0f, 1.0f);
			float X3 = Random.Range(0.0f, 2.0f*CONSTANTS.pi);
			float X4 = Random.Range(0.0f, 1.0f);
			float X5 = Random.Range(0.0f, 2.0f*CONSTANTS.pi);

			//positions
			//should I do anything to limit the number of stars that get placed very far from the center (since Plummer is infinite?)
			float psi = Mathf.Pow(Mathf.Pow(X1, -2.0f/3.0f) - 1.0f, -0.5f);
			float r = psi*rpl;
			float z = (2.0f*X2 - 1.0f)*r;
			float xy = Mathf.Sqrt(r*r - z*z);
			float x = xy*Mathf.Cos(X3);
			float y = xy*Mathf.Sin(X3);

			//velocities
			//It seems like there is a factor of 2 missing somewhere.  Using these velocities, the cluster always flies apart
			float phi = -CONSTANTS.G*Mcl/rpl*Mathf.Pow(1.0f + Mathf.Pow(r/rpl, 2.0f), -0.5f);
			vesc = Mathf.Sqrt(-2.0f*phi);
			float v = psi*vesc*vScale;
			float vz = (2.0f*X4 - 1.0f)*v;
			float vxy = Mathf.Sqrt(v*v - vz*vz);
			float vx = vxy*Mathf.Cos(X5);
			float vy = vxy*Mathf.Sin(X5);

			p.pos_data[i*3] = x + center.x;
			p.pos_data[i*3 + 1] = y + center.y;
			p.pos_data[i*3 + 2] = z + center.z;

			p.posCoM_data[i*3] = x;
			p.posCoM_data[i*3 + 1] = y;
			p.posCoM_data[i*3 + 2] = z;

			p.vel_data[i*3] = vx + vcirc.x;
			p.vel_data[i*3 + 1] = vy + vcirc.y;
			p.vel_data[i*3 + 2] = vz + vcirc.z;
			
		}

		//try this for the center of mass?
		i = Nf - 1;
		p.pos_data[i*3] = center.x;
		p.pos_data[i*3 + 1] = center.y;
		p.pos_data[i*3 + 2] = center.z;
		p.vel_data[i*3] = vcirc.x;
		p.vel_data[i*3 + 1] = vcirc.y;
		p.vel_data[i*3 + 2] = vcirc.z;
		p.mass_data[i*3] = Mcl; 
		p.rad_data[i] = 0.0f; 

		return p;
	}
	void Start () {

		cameraMover = GameObject.Find("CameraCenter").GetComponent<CameraMover>();
		bc = GameObject.Find("ButtonController").GetComponent<ButtonController>();

		//initialize all the buffers
		//something is strange about the amoutn of size needed here;
		// with floatSize = 12, this works for 10k, 15k and 30k particles, but not 20k
		pos_buf = new ComputeBuffer(NumBodies, 3*floatSize, ComputeBufferType.Default);
		posCoM_buf = new ComputeBuffer(NumBodies, 3*floatSize, ComputeBufferType.Default);
		vel_buf = new ComputeBuffer(NumBodies, 3*floatSize, ComputeBufferType.Default);
		mass_buf = new ComputeBuffer(NumBodies, 3*floatSize, ComputeBufferType.Default); //mnow, m0, mcore (currently don't have mcore)
		rad_buf = new ComputeBuffer(NumBodies, floatSize, ComputeBufferType.Default);
		lum_buf = new ComputeBuffer(NumBodies, floatSize, ComputeBufferType.Default);
		tgb_buf = new ComputeBuffer(NumBodies, 4*floatSize, ComputeBufferType.Default);
		lgb_buf = new ComputeBuffer(NumBodies, 4*floatSize, ComputeBufferType.Default);

		// These global buffers apply to every shader with these buffers defined.
		Shader.SetGlobalBuffer(Shader.PropertyToID("position"), pos_buf);
		Shader.SetGlobalBuffer(Shader.PropertyToID("positionCoM"), posCoM_buf);
		Shader.SetGlobalBuffer(Shader.PropertyToID("velocity"), vel_buf);
		Shader.SetGlobalBuffer(Shader.PropertyToID("mass"), mass_buf);
		Shader.SetGlobalBuffer(Shader.PropertyToID("radius"), rad_buf);
		Shader.SetGlobalBuffer(Shader.PropertyToID("luminosity"), lum_buf);
		Shader.SetGlobalBuffer(Shader.PropertyToID("tgb"), tgb_buf);
		Shader.SetGlobalBuffer(Shader.PropertyToID("lgb"), lgb_buf);

		//initialize the bodies
		//Particle p = initSphere(Vector3.zero, Vector3.one*rScale, Vector3.one*vScale, 0, NumBodies);
		Particle p = initPlummer(new Vector3(Rgc, 0.0f, 0.0f), Rhm, 0, NumBodies);
		//Particle p =initSimple();

		pos_buf.SetData(p.pos_data);
		posCoM_buf.SetData(p.posCoM_data);
		vel_buf.SetData(p.vel_data);
		mass_buf.SetData(p.mass_data);
		rad_buf.SetData(p.rad_data);
		lum_buf.SetData(p.lum_data);
		tgb_buf.SetData(p.tgb_data);
		lgb_buf.SetData(p.lgb_data);

	}
	


	// FixedUpdate is called many times per frame, Update is called once per frame
	// void FixedUpdate () {
	// 	if (!bc.paused)	Time.timeScale = tScale;
	// 	dt = Time.fixedDeltaTime;
	void Update () {
		//sending units of pc, MSun, Myr, but sizes in RSun
		if (!bc.paused)	Time.timeScale = tScale;
		dt = Time.deltaTime;
		time += dt;//*tScale;
		shader.SetFloat("deltaTime", dt);
		shader.SetFloat("time", time);
		shader.SetInt("NumBodies", NumBodies);
		shader.SetFloat("Softening", Softening);
		shader.SetFloat("G", CONSTANTS.G);
		shader.SetFloat("dvelMax",dvelMax); 
		shader.SetFloat("rcolFac",rcolFac); 
		shader.SetFloat("restitution", restitution);
		shader.SetFloat("neta", neta);
		shader.SetFloat("fmGal", fmGal);

		//Debug.Log(dt+" "+Time.timeScale+" "+bc.paused);

		var numberOfGroups = Mathf.CeilToInt((float) NumBodies/BLOCKSIZE);
		shader.Dispatch(shader.FindKernel("CoMCompute"), 1, 1, 1);
		shader.Dispatch(shader.FindKernel("CSMain"), numberOfGroups, 1, 1);
	}

	void OnDestroy()
	{
		pos_buf.Dispose();
		posCoM_buf.Dispose();
		vel_buf.Dispose();
		mass_buf.Dispose();
		rad_buf.Dispose();
		lum_buf.Dispose();
		tgb_buf.Dispose();
		lgb_buf.Dispose();

	}
}
