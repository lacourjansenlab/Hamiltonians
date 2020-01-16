schema = {
    "definitions": {},
    "$schema": "http://json-schema.org/draft-07/schema#",
    "$id": "http://example.com/root.json",
    "type": "object",
    "title": "The Root Schema",
    "required": [
    "N",
    "E",
    "J",
    "MaxVib",
    "wvib",
    "S",
    "sigma",
    "units"
    ],
    "properties": {
    "N": {
      "$id": "#/properties/N",
      "type": "integer",
      "title": "The N Schema",
      "default": 10,
      "examples": [
        100
      ]
    },
    "E": {
      "$id": "#/properties/E",
      "type": "integer",
      "title": "The E Schema",
      "default": 0,
      "examples": [
        10000
      ]
    },
    "J": {
      "$id": "#/properties/J",
      "type": "integer",
      "title": "The J Schema",
      "default": 100,
      "examples": [
        10
      ]
    },
    "MaxVib": {
      "$id": "#/properties/MaxVib",
      "type": "integer",
      "title": "The Maxvib Schema",
      "default": 5,
      "examples": [
        3
      ]
    },
    "wvib": {
      "$id": "#/properties/wvib",
      "type": "integer",
      "title": "The Wvib Schema",
      "default": 1400,
      "examples": [
        1400
      ]
    },
    "S": {
      "$id": "#/properties/S",
      "type": "integer",
      "title": "The S Schema",
      "default": 1,
      "examples": [
        1
      ]
    },
    "sigma": {
      "$id": "#/properties/sigma",
      "type": "integer",
      "title": "The sigma Schema",
      "default": 0,
      "examples": [
        100
      ]
    },
    "units": {
      "$id": "#/properties/units",
      "type": "string",
      "title": "The Units Schema",
      "default": "wavenumbers",
      "examples": [
        "wavenumbers"
      ],
      "pattern": "^(.*)$"
    }
    }
}
