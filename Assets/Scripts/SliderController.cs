using UnityEngine;
using UnityEngine.UI;
using System.Collections;

public class SliderController : MonoBehaviour {
	
	public string objName = "Earth";
	public string sliderClass;

	private Slider slider;
	private float value = 1.0f;
	
	private GameObject go;

	void Awake(){
		slider = GetComponent<Slider>();
		go = GameObject.Find(objName);
	}

	void Update(){
		if (!go){
			go = GameObject.Find(objName);
			slider.value = value;
			Debug.Log("not go "+go+" "+value+" "+slider.value);
		}

		if (go && slider.value != value) {
			value = slider.value;
			var script = go.GetComponent(sliderClass);
			script.SendMessage("Receiver", value);

		}

	}

}
