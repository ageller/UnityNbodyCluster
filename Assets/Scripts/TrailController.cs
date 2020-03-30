using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class TrailController : MonoBehaviour
{
	public Camera camera;

	private TrailRenderer tr;

	void Update()
	{
		//Camera's distance to parent object
		Vector3 direction = camera.transform.position - GetComponent<Rigidbody>().position;
		float distance = direction.magnitude;

		tr = GetComponent<TrailRenderer>();
		if (distance > 50.0f){
			tr.widthMultiplier = distance/200.0f;
		} else{
			tr.widthMultiplier = 0.0f;;

		}


	}

}