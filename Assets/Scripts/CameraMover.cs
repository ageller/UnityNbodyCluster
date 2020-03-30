using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class CameraMover : MonoBehaviour {

	public float epsilon = 0.0f; //softening

	private NbodyCompute nb;

	private float Mcl;
	private float galMass;
	private Vector3 pos;

	private Rigidbody rb;

	public void MclReceiver(float val){
		Mcl = val;

		Init();
	}

	Vector3 CircularXY(){
		
		//https://stackoverflow.com/questions/14845273/initial-velocity-vector-for-circular-orbit
		Vector3 dirNorm = pos/pos.magnitude;
		float cosphi = dirNorm.x;
		float sinphi = dirNorm.z;

		float initV = Mathf.Sqrt(CONSTANTS.G*(galMass + Mcl)/pos.magnitude);
		//Debug.Log("velocity "+initV+" "+pos.magnitude+" "+galMass+" "+Mcl);
		return new Vector3(-initV*sinphi, 0.0f, initV*cosphi); 


	}

	void Attract()
	{

		//to otherBody
		Vector3 direction = rb.position; //galaxy center is fixed at 0,0,0
		float distance = direction.magnitude;

		float forceMagnitude = -CONSTANTS.G*(Mcl*galMass)/(Mathf.Pow(distance, 2) + Mathf.Pow(epsilon,2));
		Vector3 force = direction.normalized*forceMagnitude;

		rb.AddForce(force);
	}

	void Init(){
		rb = GetComponent<Rigidbody>();

		NbodyCompute nb = GameObject.Find("NbodyCompute").GetComponent<NbodyCompute>();
		galMass = nb.galMass;
		pos = new Vector3(nb.Rgc, 0.0f, 0.0f); //this needs to be the same as it is set in nb
		//Debug.Log(pos+" "+Mcl);

		GetComponent<Rigidbody>().position = pos;
		GetComponent<Rigidbody>().mass = Mcl;
		GetComponent<Rigidbody>().velocity = CircularXY();

		GetComponent<TrailRenderer>().enabled = true;

	}

	void FixedUpdate(){
	//void Update(){
		Attract();
	}



}
