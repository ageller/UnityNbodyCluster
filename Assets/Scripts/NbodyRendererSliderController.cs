using UnityEngine;
using UnityEngine.UI;
using System.Collections;

public class NbodyRendererSliderController : MonoBehaviour {
	
	public string Name;

	private NbodyRenderer rend;
	private Slider slider;

	private float value = 1.0f;

	void Awake(){
		rend = GameObject.Find("NbodyCompute").GetComponent<NbodyRenderer>();
		slider = GetComponent<Slider>();
		value = Mathf.Log10((float)rend.GetType().GetField(Name).GetValue(rend));
		slider.value = value;

	}

	void Update(){

		if (slider.value != value) {
			value = slider.value;
			rend.SendMessage(Name+"Receiver", value);

		}

	}

}
