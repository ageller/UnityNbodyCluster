using UnityEngine;
using UnityEngine.UI;
using System.Collections;

public class NbodyComputeDropdownController : MonoBehaviour {
	
	public string Name;

	private NbodyCompute compute;
	private Dropdown dropdown;

	private float value = 1.0f;
	private int index = 0;

	void Awake(){
		compute = GameObject.Find("NbodyCompute").GetComponent<NbodyCompute>();
		dropdown = GetComponent<Dropdown>();
        index = dropdown.value;
        value = float.Parse(dropdown.options[index].text);
		//Debug.Log(Name+" "+value);

	}

	void Update(){

		if (dropdown.value != index) {
	        index = dropdown.value;
	        value = float.Parse(dropdown.options[index].text);
			//Debug.Log(Name+" "+value);
			compute.SendMessage(Name+"Receiver", value);
		}

	}

}
