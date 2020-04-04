// https://wiki.unity3d.com/index.php/MouseOrbitImproved
// I added an if statement to only work when clicked
// added

using UnityEngine;
using UnityEngine.Networking;
using UnityEngine.EventSystems;
using System.Collections;

//[AddComponentMenu("Camera-Control/Mouse Orbit with zoom")]
public class MouseOrbitImproved : MonoBehaviour {
 
	public Transform target;
	public Transform panMover;

	public float distance = 5.0f;
	public float xSpeed = 2.5f;
	public float ySpeed = 2.5f;
	public float panSpeed = 0.25f;

	// public float yMinLimit = -20f;
	// public float yMaxLimit = 80f;
 
	public float distanceMin = 0.001f;
	public float distanceMax = 100000;
 
	public float friction = 0.15f;

	private Rigidbody rb;
	private Camera mainCamera;

	float x = 0.0f;
	float y = 0.0f;
	float dx = 0.0f;
	float dy = 0.0f;
	float dz = 0.0f;

	//bool isHost = false;

	// Use this for initialization
	void Start (){
		mainCamera = GameObject.Find("MainCamera").GetComponent<Camera>();
		
		Vector3 angles = transform.eulerAngles;
		x = angles.y;
		y = angles.x;
 
		rb = GetComponent<Rigidbody>();
 
		// Make the rigid body not change rotation
		if (rb != null){
			rb.freezeRotation = true;
		}

	}
 
	void LateUpdate (){

		float frictionUse = friction;
		if (Input.GetKey("space")) frictionUse = 0.0f;

		if (target){
			//x,y rotation
			var mouseBtn = Input.GetMouseButton(0);
			if (mouseBtn && !EventSystem.current.IsPointerOverGameObject()){
				dx = Input.GetAxis("Mouse X")*xSpeed;
				dy = Input.GetAxis("Mouse Y")*ySpeed;

			} else {
				dx *= (1.0f - frictionUse);
				dy *= (1.0f - frictionUse);
			}
			x += dx;
			y -= dy;
			//y = ClampAngle(y, yMinLimit, yMaxLimit);
			Quaternion rotation = Quaternion.Euler(y, x, 0);

			//distance
			if (Input.GetAxis("Mouse ScrollWheel") != 0 && !EventSystem.current.IsPointerOverGameObject ()){
				dz = Input.GetAxis("Mouse ScrollWheel");
			} else {
				dz *= (1.0f - frictionUse);
			}
			//scale the rate by the distance
			float dScale = Mathf.Pow(1.0f + distance, 0.5f);
			distance = Mathf.Clamp(distance - dz*dScale, distanceMin, distanceMax);

			RaycastHit hit;
			if (Physics.Linecast (target.position, transform.position, out hit)){
				//do something in here when another object is in front of the camera?
				//distance -=  hit.distance;
			}
			Vector3 negDistance = new Vector3(0.0f, 0.0f, -distance);
			Vector3 position = rotation*negDistance + target.position;

			//apply transformations
			mainCamera.transform.rotation = rotation;
			mainCamera.transform.position = position + panMover.position;


		}

	}
 
	public static float ClampAngle(float angle, float min, float max)
	{
		if (angle < -360F)
			angle += 360F;
		if (angle > 360F)
			angle -= 360F;
		return Mathf.Clamp(angle, min, max);
	}
}