using UnityEngine;
using UnityEngine.UI;
using System.Collections;

public class NbodyRendererSliderController : MonoBehaviour {
	
	public string Name;

	private NbodyRenderer renderer;
	private Slider slider;

	private float value = 1.0f;

	void Awake(){
		renderer = GameObject.Find("NbodyCompute").GetComponent<NbodyRenderer>();
		slider = GetComponent<Slider>();
		value = Mathf.Log10((float)renderer.GetType().GetField(Name).GetValue(renderer));
		slider.value = value;

	}

	void Update(){

		if (slider.value != value) {
			value = slider.value;
			renderer.SendMessage(Name+"Receiver", value);

		}

	}

}
