using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class CameraMover : MonoBehaviour {

	public float epsilon = 0.0f; //softening

	private NbodyCompute nb;

	private float Mcl;
	//private float galMass;
	private float Rs;
	private float rho0;
	private Vector3 pos;

	private Rigidbody rb;

	public void MclReceiver(float val){
		Mcl = val;
		//Debug.Log("Initializing camera. "+Mcl);
		Init();
	}

	//this needs to be the identical function as in the NbodyCompute script
	float MWmass(float r){

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

	Vector3 CircularXY(){
		
		//https://stackoverflow.com/questions/14845273/initial-velocity-vector-for-circular-orbit
		Vector3 dirNorm = pos/pos.magnitude;
		float cosphi = dirNorm.x;
		float sinphi = dirNorm.z;
		float galMass = MWmass(pos.magnitude);
		float initV = Mathf.Sqrt(CONSTANTS.G*(galMass + Mcl)/pos.magnitude);
		//Debug.Log("velocity "+initV+" "+pos.magnitude+" "+galMass+" "+Mcl);
		return new Vector3(-initV*sinphi, 0.0f, initV*cosphi); 


	}

	void Attract()
	{

		//to otherBody
		Vector3 direction = rb.position; //galaxy center is fixed at 0,0,0
		float distance = direction.magnitude;
		float galMass = MWmass(distance);

		float forceMagnitude = -CONSTANTS.G*(Mcl*galMass)/(Mathf.Pow(distance, 2) + Mathf.Pow(epsilon,2));
		Vector3 force = direction.normalized*forceMagnitude;

		rb.AddForce(force);
	}

	void Init(){
		rb = GetComponent<Rigidbody>();

		NbodyCompute nb = GameObject.Find("NbodyCompute").GetComponent<NbodyCompute>();
		pos = new Vector3(nb.Rgc, 0.0f, 0.0f); //this needs to be the same as it is set in nb
		//Debug.Log(pos+" "+Mcl);
		float ecc = nb.ecc;

		GetComponent<Transform>().position = pos;
		GetComponent<Rigidbody>().position = pos;
		GetComponent<Rigidbody>().mass = Mcl;
		GetComponent<Rigidbody>().velocity = CircularXY()*(1.0f - ecc); //am I using this correctly?

	}

	void FixedUpdate(){
	//void Update(){
		Attract();
	}



}
