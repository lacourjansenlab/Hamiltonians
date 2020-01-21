{
  "definitions": {},
  "$schema": "http://json-schema.org/draft-07/schema#",
  "$id": "http://example.com/root.json",
  "type": "object",
  "title": "The Root Schema",
  "required": [
    "N",
    "MaxVib",
    "E",
    "J",
    "wvib",
    "S",
    "sigma",
    "oscillator_strength",
    "units"
  ],
  "properties": {
    "N": {
      "$id": "#/properties/N",
      "type": "integer",
      "title": "The N Schema",
      "default": 0,
      "examples": [
        10
      ]
    },
    "MaxVib": {
      "$id": "#/properties/MaxVib",
      "type": "integer",
      "title": "The Maxvib Schema",
      "default": 0,
      "examples": [
        5
      ]
    },
    "E": {
      "$id": "#/properties/E",
      "type": "number",
      "title": "The E Schema",
      "default": 0.0,
      "examples": [
        0.1
      ]
    },
    "J": {
      "$id": "#/properties/J",
      "type": "number",
      "title": "The J Schema",
      "default": 0.0,
      "examples": [
        100.1
      ]
    },
    "wvib": {
      "$id": "#/properties/wvib",
      "type": "number",
      "title": "The Wvib Schema",
      "default": 0.0,
      "examples": [
        1400.1
      ]
    },
    "S": {
      "$id": "#/properties/S",
      "type": "number",
      "title": "The S Schema",
      "default": 0.0,
      "examples": [
        1.1
      ]
    },
    "sigma": {
      "$id": "#/properties/sigma",
      "type": "number",
      "title": "The Sigma Schema",
      "default": 0.0,
      "examples": [
        0.1
      ]
    },
    "oscillator_strength": {
      "$id": "#/properties/oscillator_strength",
      "type": "array",
      "title": "The Oscillator_strength Schema",
      "items": {
        "$id": "#/properties/oscillator_strength/items",
        "type": "number",
        "title": "The Items Schema",
        "default": 0.0,
        "examples": [
          1.1,
          0.1,
          0.1
        ]
      }
    },
    "units": {
      "$id": "#/properties/units",
      "type": "string",
      "title": "The Units Schema",
      "default": "",
      "examples": [
        "wavenumbers"
      ],
      "pattern": "^(.*)$"
    }
  }
}
