using UnityEngine;
using UnityEngine.UI;
using System.Collections;

public class NbodyComputeSliderController : MonoBehaviour {
	
	public string Name;

	private NbodyCompute compute;
	private Slider slider;

	private float value = 1.0f;

	void Awake(){
		compute = GameObject.Find("NbodyCompute").GetComponent<NbodyCompute>();
		slider = GetComponent<Slider>();
		//can I make this more general??
		if (Name == "NumBodies"){
			value = (int)compute.GetType().GetField(Name).GetValue(compute);
		} else if (Name == "tScale") {
			value = Mathf.Log10((float)compute.GetType().GetField(Name).GetValue(compute));
		} else {
			value = (float)compute.GetType().GetField(Name).GetValue(compute);
		}
		slider.value = value;

	}

	void Update(){

		if (slider.value != value) {
			value = slider.value;
			//Debug.Log(Name+" "+value);
			compute.SendMessage(Name+"Receiver", value);

		}

	}

}
