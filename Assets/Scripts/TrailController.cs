using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class TrailController : MonoBehaviour
{
	public Camera mainCamera;

	private TrailRenderer tr;
	private bool trailOn = false;

	void Start(){
		tr = GetComponent<TrailRenderer>();
	}

	void Update()
	{
		//trying to avoid a line starting at 8500kpc and moving to the IC Rgc
		if (Time.time > 0.1 && !trailOn){
			tr.enabled = true;
			tr.emitting = true;
			tr.time = 1000.0f;
			trailOn = true;
		}

		//Camera's distance to parent object
		Vector3 direction = mainCamera.transform.position - GetComponent<Rigidbody>().position;
		float distance = direction.magnitude;
		if (distance > 50.0f){
			tr.widthMultiplier = distance/200.0f;
		} else{
			tr.widthMultiplier = 0.0f;

		}

	}

}