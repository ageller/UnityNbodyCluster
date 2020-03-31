//idea from https://forum.unity.com/threads/orbit-pan-zoom-in-one-script.11091/
using UnityEngine;
using UnityEngine.Networking;
using UnityEngine.EventSystems;
using System.Collections;

public class panMover : MonoBehaviour {
 
	private Transform cam;

	private float panSpeed;

	void Start(){
		cam = GameObject.Find("MainCamera").GetComponent<Transform>();
		panSpeed = GameObject.Find("MainCamera").GetComponent<MouseOrbitImproved>().panSpeed;
	}

	void Update(){
		if (Input.GetMouseButton(1) && !EventSystem.current.IsPointerOverGameObject()){
			transform.rotation = cam.rotation;
			transform.Translate(Vector3.right * -Input.GetAxis("Mouse X")*panSpeed);
			transform.Translate(transform.up * -Input.GetAxis("Mouse Y")*panSpeed, Space.World);
		}
	}

}