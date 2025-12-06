![Logo_long](https://github.com/ggrecow/SQAT/assets/101704265/2800e3c7-183e-4177-b011-63fc5ed3c589)

# SQAT: a sound quality analysis toolbox for MATLAB
This is the repository of **SQAT**, an open-source **S**ound **Q**uality **A**nalysis **T**oolbox for MATLAB. It contains a collection of codes implementing key metrics for **quantitative** sound quality analysis. With **SQAT** you can conduct quick quantitative sound quality analysis on any calibrated input sound file, in Pascal units. To give a transparent indication of how close the implementations are to the original models, we provide a detailed set of verification routines. Moreover, a number of example codes and exemplary sound files are provided in order to facilitate the initial use of the algorithms.  

# Sound quality metrics available in SQAT 

<!---The following psychoacoustic-based metrics are available in the `psychoacoustic_metrics` folder :
- `Loudness_ISO532_1`: Zwicker loudness model according to ISO 532-1:2017 (see [validation](validation/Loudness_ISO532_1) and [example](examples/Loudness_ISO532_1/ex_Loudness_ISO532_1.m)).
- `Sharpness_DIN45692`: Sharpness according to DIN 45692:2009 (see [validation](validation/Sharpness_DIN45692) and [example](examples/Sharpness_DIN45692/ex_Sharpness_DIN45692.m)). 
- `Roughness_Daniel1997`: Roughness model from Daniel & Weber (see [validation](validation/Roughness_Daniel1997) and [example](examples/Roughness_Daniel1997/ex_Roughness_Daniel1997.m)).  
- `FluctuationStrength_Osses2016`: Fluctuation strength model from Osses *et al.* (see [validation](validation/FluctuationStrength_Osses2016) and [example](examples/FluctuationStrength_Osses2016/ex_FluctuationStrength_Osses2016.m)).   
- `Tonality_Aures1985`: Tonality model from Aures (see [validation](validation/Tonality_Aures1985) and [example](examples/Tonality_Aures1985/ex_Tonality_Aures1985.m)).
- `PsychoacousticAnnoyance_Zwicker1999`: psychoacoustic annoyance model from Zwicker *et al.* (see [example](examples/PsychoacousticAnnoyance_Zwicker1999/ex_PsychoacousticAnnoyance_Zwicker1999.m)).
- `PsychoacousticAnnoyance_More2010`: psychoacoustic annoyance model from More (see [example](examples/PsychoacousticAnnoyance_More2010/ex_PsychoacousticAnnoyance_More2010.m)).
- `PsychoacousticAnnoyance_Di2016`: psychoacoustic annoyance model from Di *et al.* (see [example](examples/PsychoacousticAnnoyance_Di2016/ex_PsychoacousticAnnoyance_Di2016.m)). --->

The implemented metrics available in the `psychoacoustic_metrics` folder are listed in the table below, including the release on which each metric was first introduced:

| Metric  | Model | Implementation | Validation | Example | Release |
| :---: | :---: |:---: |:---: |:---: |:---: |
| Loudness  | ISO 532-1:2017 [1]  |  [link](psychoacoustic_metrics/Loudness_ISO532_1) | [link](validation/Loudness_ISO532_1) | [link](examples/Loudness_ISO532_1) | v1.0 |
| Sharpness  |  DIN 45692:2009 [2]  |  [link](psychoacoustic_metrics/Sharpness_DIN45692) | [link](validation/Sharpness_DIN45692) | [link](examples/Sharpness_DIN45692) | v1.0 |
| Roughness  | Daniel & Weber [3]  |  [link](psychoacoustic_metrics/Roughness_Daniel1997) |  [link](validation/Roughness_Daniel1997) | [link](examples/Roughness_Daniel1997) | v1.0 |
| Fluctuation Strength  | Osses *et al.* [4]  |  [link](psychoacoustic_metrics/FluctuationStrength_Osses2016) |  [link](validation/FluctuationStrength_Osses2016) | [link](examples/FluctuationStrength_Osses2016) | v1.0 |
| Tonality  | Aures [5]  |   [link](psychoacoustic_metrics/Tonality_Aures1985) |  [link](validation/Tonality_Aures1985) | [link](examples/Tonality_Aures1985) | v1.0 |
| Psychoacoustic Annoyance  | Widmann [6] (commonly misattributed to Zwicker & Fastl [7])  |  [link](psychoacoustic_metrics/PsychoacousticAnnoyance_Widmann1992) / [link](psychoacoustic_metrics/PsychoacousticAnnoyance_Zwicker1999) | - | [link](examples/PsychoacousticAnnoyance_Widmann1992) / [link](examples/PsychoacousticAnnoyance_Zwicker1999) | v1.0 |
| Psychoacoustic Annoyance  | More [8]  |  [link](psychoacoustic_metrics/PsychoacousticAnnoyance_More2010) | - | [link](examples/PsychoacousticAnnoyance_More2010) | v1.0 |
| Psychoacoustic Annoyance  | Di *et al.* [9]  |  [link](psychoacoustic_metrics/PsychoacousticAnnoyance_Di2016) | - |  [link](examples/PsychoacousticAnnoyance_Di2016) | v1.0 |
| EPNL  | FAR Part 36 [10]  |  [link](psychoacoustic_metrics/EPNL_FAR_Part36) | [link](validation/EPNL_FAR_Part36) |  [link](examples/EPNL_FAR_Part36) | v1.1 |
| Loudness  | ECMA-418-2:2025 [11]  |  [link](psychoacoustic_metrics/Loudness_ECMA418_2) | [link](validation/Loudness_ECMA418_2) |  [link](examples/Loudness_ECMA418_2) | v1.3 |
| Roughness  | ECMA-418-2:2025 [11]  |  [link](psychoacoustic_metrics/Roughness_ECMA418_2) | [link](validation/Roughness_ECMA418_2) |  [link](examples/Roughness_ECMA418_2) | v1.3 |
| Tonality  | ECMA-418-2:2025 [11]  |  [link](psychoacoustic_metrics/Tonality_ECMA418_2) | [link](validation/Tonality_ECMA418_2) |  [link](examples/Tonality_ECMA418_2) | v1.3 |

<details>
<summary><b>References</summary>
<br> 
 
[1] International Organization for Standardization. (2017). Acoustics - Methods for calculating loudness - Part 1: Zwicker method (ISO Standard No. 532-1).

[2] Deutsches Institut für Normung. (2009). Measurement technique for the simulation of the auditory sensation of sharpness (DIN Standard No. 45692).

[3] Daniel, P., & Weber, R. (1997). Psychoacoustical Roughness: Implementation of an Optimized Model. [Acta Acustica united with Acustica](https://www.ingentaconnect.com/content/dav/aaua/1997/00000083/00000001/art00020), 83(1), 113-123.

[4] Osses, A., García, R., & Kohlrausch, A. (2016). Modelling the sensation of fluctuation strength. [Proceedings of Meetings on Acoustics](https://doi.org/10.1121/2.0000410), 28(1), 050005.  

[5] Aures, W. (1985). Berechnungsverfahren für den sensorischen Wohlklang beliebiger Schallsignale (A model for calculating the sensory euphony of various sounds). Acta Acustica united with Acustica, 59(2), 130-141.

[6] Widmann, U. (1992). Ein Modell der Psychoakustischen Lästigkeit von Schallen und seine Anwendung in der Praxis der Lärmbeurteilung (A model of the psychoacoustic annoyance of sounds and its application in noise assessment practice). Doctoral thesis, Technische Universität München.

[7] Zwicker, E., & Fastl, H. (1999). Psychoacoustics: facts and models, Second edition. Springer-Verlag. DOI: [10.1007/978-3-662-09562-1](https://doi.org/10.1007/978-3-662-09562-1)

[8] More, S. R. (2010). Aircraft noise characteristics and metrics. Doctoral thesis, Purdue University. [https://docs.lib.purdue.edu/dissertations/AAI3453255/](https://docs.lib.purdue.edu/dissertations/AAI3453255/)

[9] Di, G. Q., Chen, X. W., Song, K., Zhou, B., & Pei, C. M. (2016). Improvement of Zwicker’s psychoacoustic annoyance model aiming at tonal noises. [Applied Acoustics](https://doi.org/10.1016/j.apacoust.2015.12.006), 105, 164-170.

[10] Federal Aviation Regulations. (2005). 14 CFR Parts 36 and 91, Docket No. FAA-2003-16526; Amendment No. 36-26, 91-288. [https://www.ecfr.gov/current/title-14/appendix-Appendix%20A%20to%20Part%2036](https://www.ecfr.gov/current/title-14/appendix-Appendix%20A%20to%20Part%2036) (last viewed October 30, 2023)

[11] Ecma International. (2025). Psychoacoustic metrics for ITT equipment - Part 2 (methods for describing human perception based on the Sottek Hearing Model) (Standard No. 418-2, 4th Edition/June 2025). [https://ecma-international.org/wp-content/uploads/ECMA-418-2_4th_edition_june_2025.pdf](https://ecma-international.org/wp-content/uploads/ECMA-418-2_4th_edition_june_2025.pdf) (Last viewed 16 Nov 2025)

</details>

# Sound level meter

The following **sound pressure level** (SPL) based metrics can be calculated using the codes available in `sound_level_meter` folder (see [examples](examples/sound_level_meter)):

- SPL over time
- Frequency weightings (A, B, C, D or Z)
- Time weightings (Fast, Slow or Impulse by default, but customizable)
- Equivalent SPL
- Maximum SPL
- Sound exposure level
- SPL spectrum in 1/3 octave bands
 
# Toolbox structure
The toolbox has the following directories:
- `psychoacoustic_metrics`: this directory contains a number of algorithms implementing a specific psychoacoustic metric (see [folder](psychoacoustic_metrics)). 
- `sound_level_meter`: contains scripts to obtain sound pressure levels using different frequency weightings (A, B, C, D or Z) and time weightings (Fast, Slow or Impulse) (see [folder](sound_level_meter)). 
- `utilities`: contains some scripts that are complementary to any of the toolbox functions (see [folder](utilities)).
- `examples`: an example script is provided for each metric (see [folder](examples)).
- `sound_files`: this directory hosts reference sounds in .wav format that are used mainly by the `examples` codes (see [folder](sound_files)). 
- `validation`: this directory contains scripts used to validate each algorithm. Instructions on how to run these codes are provided in each respective folder and directly on the header of the codes (see [folder](validation)). 
- `publications`: contains scripts to reproduce figures and/or tables of publications from the toolbox authors (see [folder](publications)). 

# How to use the toolbox

1. Download or clone this repository to your local computer. One way to do that is to press the button 'Code' -> Choose 'Download ZIP' and unzip somewhere).

2. After that, you need to add the relevant folders of the toolbox to the path of your MATLAB. Open and run the `startup_SQAT.m` script to automatically perform this task. In order to avoid conflicts, the `startup_SQAT.m` needs to be used every time MATLAB is (re)started. If you just want to use the metrics, you can add manually only the relevant folders to the MATLAB path (e.g., `psychoacoustic_metrics`, `sound_level_meter`, and `utilities`). 

# How to cite this repository
If you use this toolbox in your research, we would be grateful if you help us to gain visibility by citing SQAT. This is the main citation if you need to cite the toolbox repository itself:

> Felix Greco, G., Merino-Martínez, R.,  Osses, A. & Lotinga, M. J. B. (2025). SQAT: a sound quality analysis toolbox for MATLAB. Zenodo. DOI: [10.5281/zenodo.7934709](https://doi.org/10.5281/zenodo.7934709)

> [!TIP]
> **The doi above concerns the toolbox repository itself and will always resolve to the latest release. As differences between releases may occur, it is a good practice to cite the specific SQAT version being used. You can check any changes between releases [here](https://github.com/ggrecow/SQAT/releases). If you need to cite a specific release, please consult the relevant DOI [here](https://doi.org/10.5281/zenodo.7934709). If you need to cite the current SQAT release, please refer to the "Cite this repository" feature in the "About" section of this GitHub repository.**

<!-- > **Apart from the toolbox repository itself (see above), each version released has its own doi in Zenodo. As differences between releases may occur, it is a good practice to cite the specific SQAT version being used. If you need to cite the current SQAT release, please refer to the "Cite this repository" feature in the "About" section of this GitHub repository.**-->  

# Publications about SQAT

The following paper is the main work describing SQAT and the metrics available in the first release:

> Felix Greco, G., Merino-Martínez, R., Osses, A., & Langer, S. C. (2023). SQAT: a MATLAB-based toolbox for quantitative sound quality analysis. INTER-NOISE and NOISE-CON Congress and Conference Proceedings, InterNoise23, Chiba, Japan. DOI: [10.3397/IN_2023_1075](https://doi.org/10.3397/IN_2023_1075)

> Raw data and extra scripts to reproduce all the paper figures can be found [here](publications/pub_Greco2023_Internoise). 

Additionally, here's a paper by the members of the SQAT team showing three case studies where the SQAT toolbox was used to perform all analyses:

> Osses, A., Felix Greco, G., & Merino-Martínez, R. (2023). Considerations for the perceptual evaluation of steady-state and time-varying sounds using psychoacoustic metrics. 10th Convention of the European Acoustics Association (Forum Acusticum), 11-15 September 2023, Turin, Italy. DOI: [10.61782/fa.2023.0600](https://www.doi.org/10.61782/fa.2023.0600)

> Raw data and extra scripts to reproduce all the paper figures can be found [here](https://doi.org/10.5281/zenodo.7933489).

The implementation of the psychoacoustic models from ECMA-418-2 (released in v1.3) are presented and verified in the following publication:

> Lotinga, M. J. B., Torjussen, M, & Felix Greco, G. (2025). Verified implementations of the Sottek psychoacoustic Hearing Model standardised sound quality metrics (ECMA-418-2 loudness, tonality and roughness). 11th Convention of the European Acoustics Association (Forum Acusticum), 23-26 June 2025, Málaga, Spain. [(link)](https://www.researchgate.net/publication/392904348_Verified_implementations_of_the_Sottek_psychoacoustic_Hearing_Model_standardised_sound_quality_metrics_ECMA-418-2_loudness_tonality_and_roughness) 

> Raw data and extra scripts to reproduce all the paper figures can be found [here](publications/pub_Lotinga2025_Forum_Acusticum_ECMA418_2).

# Studies using SQAT

We would be very happy to know that you find SQAT useful and have used it in your own work. In this case, please reach out so we can feature your work here. 

<!-- > 2025 **-->  
<details>
<summary><b>2025</summary>

## Journal articles

> Vourakis, M., Zea, E., Karlsson, M., Andersson, N., & Etemad, S. (2024). Installation effects on axial fans: Combined aeroacoustic and psychoacoustic perspective. [Applied Acoustics](https://doi.org/10.1016/j.apacoust.2025.110872), 240.

> Lotinga, M. J. B., Green, M. C., & Toríja, A. J. (2025). Development of psychoacoustic prediction models for short-term noise annoyance responses to unmanned aircraft systems. [The Journal of the Acoustical Society of America](https://doi.org/10.1121/10.0039056), 158 (3), 2062–2082.

> Wu, C., & Redonnet, S. (2025). A simple yet efficient data-driven model for the prediction of aircraft noise impact. [Aerospace Science and Technology](https://doi.org/10.1016/j.ast.2025.110286), 163.

> Shen, Y., Bai, Y., Liu, X., & Zang, B. (2025). Drone noise reduction using serration-finlet blade design and its psychoacoustic and social impacts. [Sustainability](https://doi.org/10.3390/su17083451), 17(8), 3451.   

> Schade, S., Merino-Martínez, R., Moreau, A., Bartels, S., & Jaron, R. (2025). Psychoacoustic evaluation of different fan designs for an urban air mobility vehicle with distributed propulsion system. [The Journal of the Acoustical Society of America](https://doi.org/10.1121/10.0036228), 157 (3), 2150–2167.

> Lotinga, M. J. B., Green, M. C., & Toríja, A. J. (2025). Human perception and response to sound from unmanned aircraft systems within ambient acoustic environments. [npj Acoustics](https://doi.org/10.1038/s44384-024-00001-6), 1:2.

## Conference publications

> Lotinga, M. J. B., Green, M. C., & Toríja, A. J. (2025). Effects of exposure to unmanned aircraft systems sound: Applying machine learning and parametric clustered-data models to human response prediction. 54th International Congress & Exhibition on Noise Control Engineering (INTER-NOISE), 24-27 August 2025, São Paulo, Brazil. [(link)](https://www.researchgate.net/publication/394919259_Effects_of_exposure_to_unmanned_aircraft_systems_sound_Applying_machine_learning_and_parametric_clustered-data_models_to_human_response_prediction)

> Ferrari, G. C., Pereira Gouveia da Silva, G., & Lima Pereira, L. T. (2025). Optimizing effective perceived noise in distributed electric propulsion with neural networks and differential propeller rotation. 54th International Congress & Exhibition on Noise Control Engineering (INTER-NOISE), 24-27 August 2025, São Paulo, Brazil. [(link)](https://www.researchgate.net/publication/394937872_Optimizing_effective_perceived_noise_in_distributed_electric_propulsion_with_neural_networks_and_differential_propeller_rotation)

> Merino-Martínez, R. & Quaroni, L. N. (2025). Human response to the noise emissions of an isolated propeller under turbulent inflow conditions. 54th International Congress & Exhibition on Noise Control Engineering (INTER-NOISE), 24-27 August 2025, Sao Paulo, Brazil. [(link)](https://www.researchgate.net/publication/395129822_Human_response_to_the_noise_emissions_of_an_isolated_propeller_under_turbulent_inflow_conditions)

> Bazilinskyy, P., Alam, M. S., & Merino-Martínez, R. (2025). Pedestrian crossing behaviour in front of electric vehicles emitting synthetic sounds: A virtual reality experiment. 54th International Congress & Exhibition on Noise Control Engineering (INTER-NOISE), 24-27 August 2025, São Paulo, Brazil. [(link)](https://www.researchgate.net/publication/392237511_Pedestrian_crossing_behaviour_in_front_of_electric_vehicles_emitting_synthetic_sounds_A_virtual_reality_experiment)

> Merino-Martínez, R. & Schade, S. (2025). Psychoacoustic analysis of the perceptual influence of rotational speed fluctuations in an urban mobility vehicle with distributed ducted fans. 54th International Congress & Exhibition on Noise Control Engineering (INTER-NOISE), 24-27 August 2025, São Paulo, Brazil. [(link)](https://www.researchgate.net/publication/395129646_Psychoacoustic_analysis_of_the_perceptual_influence_of_rotational_speed_fluctuations_in_an_urban_mobility_vehicle_with_distributed_ducted_fans)

> Deutscher, B., Stalnov. O. & Ben-Gida, H. (2025). The Effect of Acoustic Detection Constraints on Optimizing Drones’ Delivery Missions. [Proceedings of the 31st AIAA/CEAS Aeroacoustics Conference](https://arc.aiaa.org/doi/10.2514/6.2025-3414).

> Podwinska, Z., Ramos-Romero, C., Green, M. C., & Toríja, A. J. (2025). The effects of time-variant characteristics of unmanned aircraft system noise on reported annoyance. 11th Convention of the European Acoustics Association (Forum Acusticum), 23-26 June 2025, Málaga, Spain. [(link)](https://www.researchgate.net/publication/393461425_The_effects_of_time-variant_characteristics_of_unmanned_aircraft_system_noise_on_reported_annoyance)

> Ramos-Romero, C., Green, M. C., Lotinga, M. J. B., & Toríja, A. J. (2025). Integrated U-space societal acceptance assessment: energy-based and perception-based acoustic metrics. 11th Convention of the European Acoustics Association (Forum Acusticum), 23-26 June 2025, Málaga, Spain. [(link)](https://www.researchgate.net/publication/393461428_Integrated_U-space_societal_acceptance_assessment_energy-based_and_perception-based_acoustic_metrics)

> Ellis, M., Green, M. C., Lotinga, M. J. B., & Toríja, A. J. (2025). Comparison of Deep Learning and Psychoacoustic Models to Predict UAV Noise Impact in Soundscapes. 11th Convention of the European Acoustics Association (Forum Acusticum), 23-26 June 2025, Málaga, Spain. [(link)](https://www.researchgate.net/publication/393231769_Comparison_of_Deep_Learning_and_Psychoacoustic_Models_to_Predict_UAV_Noise_Impact_in_Soundscapes)

> Bazilinskyy, P., Alam, M. S., & Merino-Martínez, R. (2025). Psychoacoustic assessment of synthetic sounds for electric vehicles in a virtual reality experiment. 11th Convention of the European Acoustics Association (Forum Acusticum), 23-26 June 2025, Málaga, Spain. [(link)](https://www.researchgate.net/publication/390563218_Psychoacoustic_assessment_of_synthetic_sounds_for_electric_vehicles_in_a_virtual_reality_experiment)

> Merino-Martínez, R. & Buzeţelu, V. S. (2025). Aircraft noise-induced annoyance analysis using psychoacoustic listening experiments. 11th Convention of the European Acoustics Association (Forum Acusticum), 23-26 June 2025, Málaga, Spain. [(link)](https://www.researchgate.net/publication/393163550_Aircraft_noise-induced_annoyance_analysis_using_psychoacoustic_listening_experiments)

> Priboi, S.A. & Merino-Martínez, R. (2025). Evaluation of audio-visual parameters in the perceived aircraft noise annoyance using virtual reality experiments. 11th Convention of the European Acoustics Association (Forum Acusticum), 23-26 June 2025, Málaga, Spain. [(link)](https://www.researchgate.net/publication/393163554_Evaluation_of_audio-visual_parameters_in_the_perceived_aircraft_noise_annoyance_using_virtual_reality_experiments)

> Merino-Martínez, R. & Quaroni, L. N. (2025). Psychoacoustic characterization of an isolated propeller at different inflow turbulence conditions and collective pitch angles. 11th Convention of the European Acoustics Association (Forum Acusticum), 23-26 June 2025, Málaga, Spain. [(link)](https://www.researchgate.net/publication/393163829_Psychoacoustic_characterization_of_an_isolated_propeller_at_different_inflow_turbulence_conditions_and_collective_pitch_angles)

> Lladó, P., Neidhardt, A., Brinkmann, F., & de Sena, E. (2025). Spatial audio models' inventory to cover the attributes from the spatial audio quality inventory. 11th Convention of the European Acoustics Association (Forum Acusticum), 23-26 June 2025, Málaga, Spain. [(link)](https://www.researchgate.net/profile/Pedro-Llado/publication/393003782_Spatial_audio_models'_inventory_to_cover_the_attributes_from_the_spatial_audio_quality_inventory/links/685bd45799d2ce32c1cac97f/Spatial-audio-models-inventory-to-cover-the-attributes-from-the-spatial-audio-quality-inventory.pdf) 

> Pockelé, J. S. & Merino-Martínez, R. (2025). Perceived Noise Impact of Transitioning Towards Larger Wind Turbines Using Auralisations. 11th Ed. International Conferences on Wind Turbine Noise, Copenhagen, Denmark. [(link)](https://www.researchgate.net/publication/393400798_Perceived_Noise_Impact_of_Transitioning_Towards_Larger_Wind_Turbines_Using_Auralisations)

> Pockelé, J. S. & Merino-Martínez, R. (2025). Influence of Ambient Noise in Sound Quality Assessment of Auralised Wind Turbine Noise. 11th Ed. International Conferences on Wind Turbine Noise, Copenhagen, Denmark. [(link)](https://www.researchgate.net/publication/393401282_Influence_of_Ambient_Noise_in_Sound_Quality_Assessment_of_Auralised_Wind_Turbine_Noise)

<br> 
</details>

<!-- > 2024 **-->  

<details>
<summary>2024</summary>

 ## Journal articles

> Schmidt, H., Yupa-Villanueva, R. M., Ragni, D., Merino-Martínez, R., van Gool, P., & Schmehl, R. (2024). Exploring noise annoyance and sound quality for airborne wind energy systems: insights from a listening experiment. [Wind Energy Science](https://doi.org/10.5194/wes-10-579-2025), 10, 579–595.

> Kawai, C., Jäggi, J., Georgiou, F., Meister, J., Pieren, R., & Schäffer, B. (2024). Short-term noise annoyance towards drones and other transportation noise sources: A laboratory study. [The Journal of the Acoustical Society of America](https://doi.org/10.1121/10.0032386), 156 (4), 2578–2595.

> Louwers, G., Pont, S., Gommers, D., van der Heide, E., & Özcan, E. (2024). Sonic ambiances through fundamental needs: An approach on soundscape interventions for intensive care patients, [The Journal of the Acoustical Society of America](https://doi.org/10.1121/10.0030470), 156 (4), 2376–2394.

> Brandetti, L., Mulders, S. P., Merino-Martínez, R., Watson, S., & van Wingerden, J.-W. (2024). Multi-objective calibration of vertical-axis wind turbine controllers: balancing aero-servo-elastic performance and noise. [Wind Energy Science](https://doi.org/10.5194/wes-9-471-2024), 9, 471-493.

## Conference publications

> Lotinga, M. J. B., Green, M. C., & Toríja, A. J. (2024). How do flight operations and ambient acoustic environments influence noticeability and noise annoyance associated with unmanned aircraft systems? [Quiet Drones 2024 conference](https://www.researchgate.net/publication/383915149_How_do_flight_operations_and_ambient_acoustic_environments_influence_noticeability_and_noise_annoyance_associated_with_unmanned_aircraft_systems).  

> Merino-Martínez, R., Yupa Villanueva, R. M., von den Hoff, B., & Pockelé, J. S. (2024). Human response to the flyover noise of different types of drones recorded in field measurements. [Quiet Drones 2024 conference](https://www.researchgate.net/publication/384065422_Human_response_to_the_flyover_noise_of_different_types_of_drones_recorded_in_field_measurements).

> Snellen, M., Merino-Martinez, R., Altena, A., Amiri-Simkooei, A., Andino Cappagli, C.I., Morin, A., Quaroni, L.N., Yunus, F. & Yupa-Villanueva, R.M. (2025) Research on drone and urban air mobility noise: Measurement, modelling, and human perception. [Quiet Drones 2024 conference](https://www.researchgate.net/publication/386098572_Research_on_drone_and_urban_air_mobility_noise_Measurement_modelling_and_human_perception_Session).

> Georgiou, F., Schäffer, B., Heusser, A., & Pieren, R. (2024). Prediction of Noise Annoyance of Air Vehicle Flyovers Using Psychoacoustic Models. [Proceedings of the 30th International Congress on Sound and Vibration (ICSV)](https://www.researchgate.net/publication/384355002_PREDICTION_OF_NOISE_ANNOYANCE_OF_AIR_VEHICLE_FLY-_OVERS_USING_PSYCHOACOUSTIC_MODELS).

> Yupa Villanueva, R. M., Merino-Martínez, R., Andino Cappagli, C. I., Altena, A., & Snellen, M. (2024). Effect of Unmanned Aerial Vehicle Configurations on the Acoustic and Psychoacoustic Signatures. [Proceedings of the 30th International Congress on Sound and Vibration (ICSV)](https://www.researchgate.net/publication/385660330_EFFECT_OF_UNMANNED_AERIAL_VEHICLE_CONFIGURATIONS_ON_THE_ACOUSTIC_AND_PSYCHOACOUSTIC_SIGNATURES).

> von den Hoff, B., Merino-Martínez, R., Yupa Villanueva, R. M., & Snellen, M. (2024). Noise Emissions and Noise Annoyance of a Single-Propeller Electric Aircraft During Flyover. [Proceedings of the 30th International Congress on Sound and Vibration (ICSV)](https://www.researchgate.net/publication/382695675_NOISE_EMISSIONS_AND_NOISE_ANNOYANCE_OF_A_SINGLE-PROPELLER_ELECTRIC_AIRCRAFT_DURING_FLYOVER).

> Pockelé, J. S. & Merino-Martínez, R. (2024). Psychoacoustic Evaluation of Modelled Wind Turbine Noise. [Proceedings of the 30th International Congress on Sound and Vibration (ICSV)](https://www.researchgate.net/publication/382256126_PSYCHOACOUSTIC_EVALUATION_OF_MODELLED_WIND_TURBINE_NOISE).

> Merino-Martínez, R., Ben-Gida, H., & Snellen, M. (2024). Psychoacoustic Evaluation of an Optimized Low-Noise Drone Propeller Design. [Proceedings of the 30th International Congress on Sound and Vibration (ICSV)](https://www.researchgate.net/publication/382255303_PSYCHOACOUSTIC_EVALUATION_OF_AN_OPTIMIZED_LOW-NOISE_DRONE_PROPELLER_DESIGN).
 
> Yupa Villanueva, R.M., Merino-Martínez, R., Altena, A., & Snellen, M. (2024). Psychoacoustic Characterization of Multirotor Drones in Realistic Flyover Maneuvers. [Proceedings of the 30th AIAA/CEAS Aeroacoustics Conference](https://arc.aiaa.org/doi/10.2514/6.2024-3015).

> Thoma, E.M., Merino-Martínez, R., Grönstedt, T., & Zhao, X. (2024). Noise From Flight Procedure Designed With Statistical Wind: Auralization and Psychoacoustic Evaluation. [Proceedings of the 30th AIAA/CEAS Aeroacoustics Conference](https://arc.aiaa.org/doi/10.2514/6.2024-3017).

> Schade, S., Merino-Martínez, R., Ratei, P., Bartels, S., Jaron, R., & Moreau, A. (2024). Initial Study on the Impact of Speed Fluctuations on the Psychoacoustic Characteristics of a Distributed Propulsion System with Ducted Fans. [Proceedings of the 30th AIAA/CEAS Aeroacoustics Conference](https://arc.aiaa.org/doi/10.2514/6.2024-3273).

> Monteiro, F.d. N., Merino-Martínez, R., & Lima Pereira, L.T. (2024). Psychoacoustic Evaluation of an Array of Distributed Propellers Under Synchrophasing Operation. [Proceedings of the 30th AIAA/CEAS Aeroacoustics Conference](https://arc.aiaa.org/doi/10.2514/6.2024-3321).

> Merino-Martínez, R., Besnea, I., von den Hoff, B., & Snellen, M. (2024). Psychoacoustic Analysis of the Noise Emissions from the Airbus A320 Aircraft Family and its Nose Landing Gear System. [Proceedings of the 30th AIAA/CEAS Aeroacoustics Conference](https://arc.aiaa.org/doi/10.2514/6.2024-3398).

> Knuth, D., Ring, T. P., & Langer, S. C. (2024). Comparing auralizations and measurements of vibrating plates with physical and psychoacoustic metrics. [Proceedings of 50. Jahrestagung für Akustik (DAGA)](https://pub.dega-akustik.de/DAGA_2024/files/upload/paper/176.pdf).
<br>

</details>
 
<!-- > 2023 **-->  
<details>
<summary>2023</summary>

## Conference publications
 
> Knuth, D., Ring, T. P., & Langer, S. C. (2023). Utilizing auralization to investigate psychoacoustic perception of vibrating structures. [Proceedings of 49. Jahrestagung für Akustik (DAGA)](https://pub.dega-akustik.de/DAGA_2023/data/articles/000414.pdf).
<br> 
</details>

<!-- > thesis **-->  
<details>
<summary>Thesis</summary>
<br>

> Schade, S. (2025). Design of low-noise fan engines for urban air mobility and sound quality analysis using virtual flyovers. Doctoral thesis, Technische Universität Berlin. DOI: [10.14279/depositonce-24074](https://doi.org/10.14279/depositonce-24074)

> Gan, Z. F. (2025). Time-Varying Noise of Electric Multirotor Aircraft. Doctoral thesis, Pennsylvania State University. [(link)](https://etda.libraries.psu.edu/catalog/25752zug117) 

> Priboi, S. A. (2025). Evaluation of audio-visual parameters in the perceived aircraft noise annoyance using virtual reality experiments. Master thesis, Delft University of Technology. [(link)](https://resolver.tudelft.nl/uuid:70a31852-ddaa-4e5c-90d0-d9c47352da6c)

> Buzeţelu, V. S. (2025). Aircraft-Induced Psychoacoustic Annoyance Quantification Using Artificial Intelligence. Master thesis, Delft University of Technology. [(link)](https://resolver.tudelft.nl/uuid:fa3ddab2-3a9b-4ad6-97d1-c3ce0a4d9678)

> Taniguchi, R. (2024). 音質評価指標を用いた感覚的快さと感覚的快くなさの評価に関する調査 (Survey on the evaluation of sensory pleasantness and unpleasantness using sound quality evaluation indexes). Master thesis, Japan Advanced Institute of Science and Technology. [(link)](http://hdl.handle.net/10119/18909)

> Pockelé, J. S. (2023). Auralisation of modelled wind turbine noise for psychoacoustic listening experiments: development and validation of the wind turbine auralisation tool WinTAur. Master thesis, Delft University of Technology. [(link)](http://resolver.tudelft.nl/uuid:cc9e67b4-6bde-4114-97c0-43b11b4a48ef)

</details>

# Toolbox history

Motivated by the limited access to well-documented and validated open-source implementations of perceptually-inspired metrics, [Gil Felix Greco](https://www.linkedin.com/in/gil-felix-greco-363985101/) initiated in the end of 2019 a compilation of algorithms to be used in the evaluation of environmental aircraft noise. In 2021, this set of algorithms was compiled into a MATLAB toolbox, which was named for the first time as a sound quality analysis toolbox (SQAT), and applied to analyze the sound quality of a novel aircraft concept (see publication below). This preliminary SQAT version (not hosted in this repository) required a very specific format of input parameters provided by an aircraft noise prediction research software (PANAM/DLR), suited to the needs of the aforementioned publication.

> Felix Greco, G., Bertsch, L., Ring, T. P., & Langer, S. C. (2021). Sound quality assessment of a medium-range aircraft with enhanced fan-noise shielding design. [CEAS Aeronautical Journal](https://doi.org/10.1007/s13272-021-00515-9) 12, 481–493.

The actual conception of the SQAT toolbox as published in this repository occurred in October 2022, after inspiring discussions during the ICA conference in South Korea, where Gil met Alejandro and Roberto. [Alejandro Osses](https://www.linkedin.com/in/alejandro-osses-10039883/) has been actively contributing to the hearing research community since 2014 and as a developer of the AMT toolbox since 2020 (amtoolbox.org). He provided know-how about code structure and the use of GitHub. [Roberto Merino-Martínez](https://www.linkedin.com/in/roberto-merino-martinez/) has been active in the use and development of sound quality metrics and their actual validation through listening experiments, also contributing some of his code implementations to the team. As a result of Gil's previous efforts and the newly joined forces, we were able to publish the first release of the toolbox, SQAT v1.0 in May 2023. In 2025, [Mike Lotinga](https://www.linkedin.com/in/michael-lotinga-6890328/) joined the team of maintainers after contributing to SQAT v1.3 with the implementation of ECMA-418-2 models.


Within the SQAT team, we are committed to carefully maintaining the toolbox in the long term, thus supporting its use as a reliable and readily accessible tool. The inclusion of new metrics in the future is foreseen as long as enough verification and documentation related to their implementations are (or can be made) available. With this strong requirement, we will only publish implementations that are mature enough for general use, without necessarily requiring further updates after a code is released. We are committed to ensuring backward compatibility of the codes across SQAT releases. However, small unforeseen issues will be fixed as quickly as possible. If that is the case, any changes will be clearly stated it in the documentation. For this reason, we recommend that users always inform which version of SQAT has been used, and to use the most up-to-date release of the toolbox.

# Contact

If you would like to get in touch to report a bug, make suggestions or ask a question, please contact us on GitHub by opening an issue.

# Licensing
<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0 International License</a>.